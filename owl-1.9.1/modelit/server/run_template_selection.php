<?php
include_once("inc/settings.php");
session_start();

function create_output_dir($job_name) {
	$j = escapeshellcmd($job_name);
	$output_dir = get_output_dir($job_name);
	$cmd = "mkdir ".$output_dir;
	exec($cmd);
}

# writes two files
# 1: $job_name.fa - subsequence of $sequence from $seq_from to $seq_to
# 2: $job_name.full.fa - full $sequence with above subsequence in upper and rest in lower case
# fasta header of 1 is the job name
# fasta header of 2 is the subsequence index in a format like: >123-456
# assumes that $sequence contains only amino acid characters and $seq_from and $seq_to are valid indices counting from 1
function create_seq_files($job_name, $fullseq, $subseq, $seq_from, $seq_to) {
	$seq_file = get_seq_file($job_name);
	$full_seq_file = get_full_seq_file($job_name);
	
	// write subsequence
	$fh = fopen($seq_file, 'w') or die("can't open file ".$seq_file);
	fwrite($fh, ">".$job_name."\n");
	fwrite($fh, $subseq."\n");
	fclose($fh);
	
	// write full sequence
	$fh = fopen($full_seq_file, 'w') or die("can't open file ".$full_seq_file);
	fwrite($fh, ">".$seq_from.'-'.$seq_to."\n");
	fwrite($fh, $fullseq."\n");
	fclose($fh);
}

function get_seq_file($job_name) {
	$j = escapeshellcmd($job_name);
	return $GLOBALS['server_root']."results/".$j."/".$j.".fa";		
}

function get_full_seq_file($job_name) {
	$j = escapeshellcmd($job_name);
	return $GLOBALS['server_root']."results/".$j."/".$j.".full.fa";		
}

function get_error_log_url($job_name) {
	$j = escapeshellcmd($job_name);
	return "results/".$j."/".$j.".err.log";		
}

function get_error_log($job_name) {
	$j = escapeshellcmd($job_name);
	return $GLOBALS['server_root']."results/".$j."/".$j.".err.log";	
}

function get_out_log($job_name) {
	$j = escapeshellcmd($job_name);
	return $GLOBALS['server_root']."results/".$j."/".$j.".out.log";	
}

function get_ts_job_name($job_name) {
	return 'TS-'.$job_name;
}

function get_ts_lock_file($job_name) {
	return get_result_file($job_name, ".TS.lock");
}

// writes a lock file indicating that template selection was started
function write_ts_lock_file($job_name) {
	touch(get_ts_lock_file($job_name));
}

function run_template_selection($job_name) {
	global $use_parallel_blast;
	global $num_blast_cpus;
    global $num_templates_in_report_file;
    global $eval_cutoff_template_selection;
    global $num_psiblast_iterations;
	write_ts_lock_file($job_name);
	$output_dir = get_output_dir($job_name);
	$seq_file = get_seq_file($job_name);
	$error_log = get_error_log($job_name);
	$output_log = get_out_log($job_name);
	
	$n = $num_blast_cpus;
	if($use_parallel_blast) {
		# multi threaded blast search (default: 8 cpus)
		$cmd = "qsub -pe threaded $n -V -N ".get_ts_job_name($job_name)." -e ".get_error_log($job_name)." -o ".get_out_log($job_name)." -q all.q /project/StruPPi/CASP8/scripts/templateSelection -i $seq_file -a $n -o $output_dir -j $num_psiblast_iterations -e $eval_cutoff_template_selection -x $num_templates_in_report_file";
	} else {
		# single threaded blast search
		$cmd = "qsub -V -N $job_name -e ".get_output_dir($job_name)." -o ".get_output_dir($job_name)." -q all.q /project/StruPPi/CASP8/scripts/templateSelection -i $seq_file -o $output_dir -j $num_psiblast_iterations -e $eval_cutoff_template_selection -x $num_templates_in_report_file";
	}
	
	$retval = 0;
	system($cmd,$retval);
	echo '<p></p>';
	#echo '<p>'.$cmd.'</p>';
	return $retval;
}

if(!isset($_SESSION['job_name']) || !isset($_SESSION['sequence'])) {
	header("Location: templates.php");
} else {
	$job_name = $_SESSION['job_name'];
	$sequence = $_SESSION['sequence'];	// full sequence in upper case
	$seq_from = 1;
	$seq_to = strlen($sequence);
	if(array_key_exists('seqfrom',$_SESSION)) $seq_from = $_SESSION['seqfrom'];
	if(array_key_exists('seqto',$_SESSION)) $seq_to = $_SESSION['seqto'];
	
	$subseq = substr($sequence, $seq_from - 1, $seq_to-$seq_from + 1);
	$fullseq = strtolower(substr($sequence,0,$seq_from-1)).strtoupper($subseq).strtolower(substr($sequence,$seq_to));
}
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
  <title>Structure prediction server</title>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
  <link rel="stylesheet" type="text/css" href="css/default.css" />
</head>
<body>
<?php include_once("inc/owl_header.inc.php") ?>
<h3>Template search</h3>
<?php
	
	echo "<p><b>Job name:</b> ".$job_name."</p>";
	
	echo '<p><b>Target sequence: </b>'.$subseq.'</p>';	
	
	echo '<p><b>Position in full sequence: </b>'.$seq_from.'-'.$seq_to.'</p>';	
	
	// create output dir
	create_output_dir($job_name);
	
	// write sequence to output dir
	create_seq_files($job_name, $fullseq, $subseq, $seq_from, $seq_to);
	
	if(file_exists(get_ts_lock_file($job_name))) {
		echo '<p>Job already submitted.</p>';
	} else {
		// run templateSelection with output dir and sequence file
		$result = run_template_selection($job_name);
		if( $result != 0) {
			echo '<p><font color="red">An error occured </font>(exit status='.$result.'). See the <a href="'.get_error_log_url($job_name).'">error log</a>.</p>';		
		}	
	}
?>

<div class="buttons">
<a href="wait_for_templates.php?job_name=<? echo $job_name; ?>"><img src="images/eye.png" alt="" /> Check results</a>
</div>

</body>
</html>
