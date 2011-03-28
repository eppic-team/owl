<?php
include_once("inc/settings.php");

if(!isset($_GET['job_name'])) {
	error('No job name specified.');
} else {
	$job_name = $_GET['job_name'];
}

function get_result_dir($job_name) {
	$j = escapeshellcmd($job_name);
	return $GLOBALS['server_root']."results/".$j."/";	
}

function get_template_file($job_name) {
	return get_result_file($job_name, ".man.templates");
}

function get_seq_file($job_name) {
	return get_result_file($job_name, ".fa");
}

function get_mb_job_name($job_name) {
	return 'MB-'.$job_name;
}

function run_modeling($job_name) {
	global $use_parallel_tinker;
	write_mb_lock_file($job_name);
	$t_file = get_template_file($job_name);
	$s_file = get_seq_file($job_name);
	$o_dir = get_result_dir($job_name);
	
	# old shell-script
	#$cmd = sprintf("qsub -V -q all.q -N %s -o %s -e %s /project/StruPPi/CASP8/scripts/runall.sh -i %s -o %s -t %s",get_mb_job_name($job_name), $o_dir, $o_dir, $s_file, $o_dir, $t_file);
		
	if($use_parallel_tinker) {
	    # brand new parallel tinker java code
		$cmd = sprintf("qsub -V -q all.q -N %s -o %s -e %s /project/StruPPi/CASP8/scripts/model_it -A -i %s -o %s -t %s",get_mb_job_name($job_name), $o_dir, $o_dir, $s_file, $o_dir, $t_file);			
	} else {
		# standard java code		
		$cmd = sprintf("qsub -V -q all.q -N %s -o %s -e %s /project/StruPPi/CASP8/scripts/model_it -i %s -o %s -t %s",get_mb_job_name($job_name), $o_dir, $o_dir, $s_file, $o_dir, $t_file);
	}
	
	$retval = 0;
	system($cmd,$retval);
	#echo($cmd);
	echo '<p></p>';
	return $retval;
}

function get_mb_lock_file($job_name) {
	return get_result_file($job_name, ".MB.lock");
}

// writes a lock file indicating that modeling was started so that it won't be resubmitted when the browser's
// 'back' button is used. The creation time of the file indicates when the job was started. The finishing time
// can be taken from the creation time of the pdb file.
function write_mb_lock_file($job_name) {
	touch(get_mb_lock_file($job_name));
}

# make sure that template file exists, otherwise can't model
if(!file_exists(get_template_file($job_name))) {
	error('Template file for job'.$job_name.' not found.');
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
<h3>Run structure modelling</h3>
<?php
	
	echo "<p><b>Job name:</b> ".$job_name."</p>";
	
	echo '<p><b>Templates:</b><br/>';
	$templates = file(get_template_file($job_name));
	foreach($templates as $key => $val) {
		echo $val."<br/>";
	}
	echo '</p>';
		
	// create output dir
	#echo "<p>Creating output directory...<br/>";
	#create_output_dir($job_name);
	
	if(file_exists(get_mb_lock_file($job_name))) {
		echo '<p>Job already submitted.</p>';
	} else {
		
		// run templateSelection with output dir and sequence file
		$result = run_modeling($job_name);
		if( $result != 0) {
			echo '<p><font color="red">An error occured </font>(exit status='.$result.'). See the <a href="'.get_error_log_url($job_name).'">error log</a>.</p>';		
		}		
	}


?>

<div class="buttons">
<a href="wait_for_modeling_result.php?job_name=<? echo $job_name; ?>"><img src="images/eye.png" alt="" /> Check results</a>
</div>

</body>
</html>