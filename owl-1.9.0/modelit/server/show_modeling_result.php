<?php
include_once("inc/settings.php");

if(isset($_GET['job_name'])) {
	$job_name = $_GET['job_name'];
}

if(!isset($job_name)) {
	$errors['no_job_name'] = 1;
}

function get_output_url($job_name) {
	return "results/".$job_name."/";
}

function get_orig_pdb_file($job_name) {
	return get_result_file($job_name, ".reconstructed.pdb");
}

function get_renum_pdb_file($job_name) {
	return get_result_file($job_name, ".renum.pdb");
}

function get_pdb_url($job_name) {
	if(file_exists(get_renum_pdb_file($job_name))) {
		return get_output_url($job_name).$job_name.".renum.pdb";
	} else {
		return get_output_url($job_name).$job_name.".reconstructed.pdb";
	}
}

function get_img_file($job_name) {
	return get_result_file($job_name, ".png");
}

function get_img_url($job_name) {
	return get_output_url($job_name).$job_name.".png";
}

function get_blast_img($job_name) {
	return get_output_url($job_name)."blast.png";	
}

function get_psiblast_img($job_name) {
	return get_output_url($job_name)."psiblast.png";		
}

function get_ss_compare_file($job_name) {
	return get_result_file($job_name, ".ss.compare");
}

function get_ss_compare_img_file($job_name) {
	return get_result_file($job_name, ".ss.compare.png");
}

function get_ss_compare_img_url($job_name) {
	return get_output_url($job_name).$job_name.".ss.compare.png";
}

function get_ss_pred_file($job_name) {
	return get_result_file($job_name, ".horiz");
}

function get_template_file($job_name) {
	return get_result_file($job_name, ".man.templates");
}

function read_templates($job_name) {
	$result = array();
	if(file_exists(get_template_file($job_name))) {
		$result = file(get_template_file($job_name), FILE_SKIP_EMPTY_LINES);
	}
	return $result;
}

function get_seq_file($job_name) {
	return get_result_file($job_name, ".fa");
}

function get_full_seq_file($job_name) {
	return get_result_file($job_name, ".full.fa");
}

function read_sequence($job_name) {
	if(file_exists(get_seq_file($job_name))) {
		$lines = file(get_seq_file($job_name));
		if(count($lines) < 2) {
			return "";
		} else {
			$seq = $lines[1];
			return $seq;
		}
	} else {
		return "";
	}
}

function read_full_sequence($job_name) {
	if(file_exists(get_full_seq_file($job_name))) {
		$lines = file(get_full_seq_file($job_name));
		if(count($lines) < 2) {
			return "";
		} else {
			$seq = $lines[1];
			return $seq;
		}
	} else {
		return "";
	}
}

function read_subseq_interval($job_name) {
	if(file_exists(get_full_seq_file($job_name))) {
		$lines = file(get_full_seq_file($job_name));
		if(count($lines) < 2) {
			return "";
		} else {
			$header = $lines[0];
			return trim(substr($header,1));
		}
	} else {
		return "";
	}
}

function print_seq($seq) {
	echo '<div class="sequence">';
	$upper = 1;
	for($i = 0; $i < strlen($seq); $i++) {
		if($i > 0 && $i % 100 == 0) {
			echo '<br/>';
		}
		if($upper && strtoupper($seq[$i])!=$seq[$i]) {
			$upper = 0;
			echo '<span class="fadeseq">';
		}
		if(!$upper && strtoupper($seq[$i])==$seq[$i]) {
			$upper = 1;
			echo '</span>';
		}
		echo $seq[$i];
	}
	echo '</div>';
}

function create_png_preview($job_name) {
	$j = escapeshellcmd($job_name);
	$pdb_file = get_orig_pdb_file($job_name);
	$png_file = get_img_file($job_name);
	
	$width = 480;
	$height = 360;
	$cmd = $GLOBALS['script_dir']."pdb2png.sh ".$pdb_file." ".$width." ".$height." ".$png_file;
	exec($cmd);	
}

function create_ss_compare_img($job_name) {
	$cmp_file = escapeshellcmd(get_ss_compare_file($job_name));
	$img_file = escapeshellcmd(get_ss_compare_img_file($job_name));
	$width = 880;
	$cmd = $GLOBALS['script_dir']."render_ss_compare.pl $cmp_file $width > $img_file";
	$result = exec($cmd);
	//echo $cmd."<br/>";
}

function renum_pdb($job_name, $offset) {
	$offs = $offset + 1; # needed because script works unintuitively, taking new starting residue instead of offset
	$old_file = escapeshellcmd(get_orig_pdb_file($job_name));
	$new_file = escapeshellcmd(get_renum_pdb_file($job_name));
	$chain = 'A';
	$cmd = $GLOBALS['script_dir']."renumpdb.py ".$offs." ".$old_file." ".$new_file." ".$chain;
	$result = exec($cmd);
}

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
  <title>Structure prediction server</title>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
  <link rel="stylesheet" type="text/css" href="css/default.css" />

  <script type="text/javascript">
      
    function toggle_fieldset(infobox_id) {
    	var ib_obj = document.getElementById(infobox_id);
    	var fs_obj = ib_obj.getElementsByTagName('fieldset')[0];
    	var div_obj = fs_obj.getElementsByTagName('div')[0];
    	var d = div_obj.style.display;
    	if(d == "none") {
    		div_obj.style.display = "";
    		fs_obj.style.border = ''; //"1px solid #CCCCCC";
    		fs_obj.style.margin = ''; //"top:5px bottom:5px"
    	} else {
    		div_obj.style.display = "none";
    		fs_obj.style.borderBottom = "none";
    		fs_obj.style.borderLeft = "none";
    		fs_obj.style.borderRight = "none";    		
    		fs_obj.style.marginBottom = "-10px";
       	}
    }
  
  </script>

</head>

<body>
<?php include_once("inc/owl_header.inc.php") ?>

<h3>Modeling result</h3>

<? # job name ?>
<div class="infobox" id="ib_name">
<form>
<fieldset>
<legend  onclick="javascript:toggle_fieldset('ib_name')">Jobname</legend>
<div>
<? echo $job_name; ?>
</div>
</fieldset>
</form>
</div>

<? # full sequence ?>
<? $full_seq = read_full_sequence($job_name); 
if(strlen($full_seq) > 0) { ?>
<div class="infobox" id="ib_full_seq">
<form><fieldset>
<legend onclick="javascript:toggle_fieldset('ib_full_seq')" >Full sequence</legend>
<div>
<? print_seq($full_seq); ?>
</div>
</fieldset></form>
</div>
<? } ?>

<? # sequence ?>
<div class="infobox" id="ib_seq">
<form>
<fieldset>
<?  $intv = read_subseq_interval($job_name); if(strlen($intv) > 0) { $intv = '('.$intv.')'; } ?>
<legend onclick="javascript:toggle_fieldset('ib_seq')" >Target sequence <? echo $intv; ?> </legend>
<div>
<? print_seq(read_sequence($job_name)); ?>
</div>
</fieldset>
</form>
</div>

<? # templates ?>
<div class="infobox" id="ib_temps">
<fieldset>
<legend onclick="javascript:toggle_fieldset('ib_temps')">Templates used for modeling</legend>
<div>
<?php
$temps = read_templates($job_name);
if(count($temps) > 0) {
	foreach($temps as $temp) {
		$code = rtrim($temp);
		$pdb = substr($code, 0, 4);
		echo '<a href="http://pdbwiki.org/index.php/'.$pdb.'">'.rtrim($temp).'</a> ';
	}
} else {
	echo 'Template file not found';
}
?>
</div>
</fieldset>
</div>

<? # predicted structure ?>
<?php  
	$intv = read_subseq_interval($job_name);
	if(strlen($intv) > 0) {
		$start_end = split("-",$intv);
		$offset = $start_end[0] - 1;
		if($offset > 0) {
			renum_pdb($job_name, $offset);
		}
	} 
?>
<div class="infobox" id="ib_pdb">
<form>
<fieldset>
<legend onclick="javascript:toggle_fieldset('ib_pdb')">Predicted structure</legend>
<div>
<?php create_png_preview($job_name); ?>
<A href="<?php echo get_pdb_url($job_name); ?>"><IMG src="<? echo get_img_url($job_name); ?>" border="0"></A>
</div>
</fieldset>
</form>
</div>

<? # secondary structure 
if($GLOBALS['show_ss_on_result_page']) {
?>
<div class="infobox" id="ib_ss">
<form>
<fieldset>
<legend onclick="javascript:toggle_fieldset('ib_ss')">Secondary structure</legend>
<div>
<?
	create_ss_compare_img($job_name);
	if(file_exists(get_ss_compare_file($job_name))) {
		echo '<img src="'.get_ss_compare_img_url($job_name).'" alt="">';
	} else {
//	if(file_exists(get_ss_pred_file($job_name))) {
//		$cont = file_get_contents(get_ss_pred_file($job_name));
//		echo '<pre>'.$cont.'</pre>';				
//	} else {
		echo "Secondary structure file not found!";
	}
?>
</div>
</fieldset>
</form>
</div>
<? } ?>

<div class="buttons">
<p><a href="<?php echo get_pdb_url($job_name); ?>"><img src="images/pdb.png" alt="" /> Download pdb file</a></p>
<a href="results.php"><img src="images/application_view_list.png" alt="" /> View result list</a>
<p><a href="index.php"><img src="images/new.png" alt="" /> Start new job</A></p>
</div>
<div style="clear:both;"></div>
<p></p>

</body>
</html>