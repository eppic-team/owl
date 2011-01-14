<?php
include_once("inc/settings.php");

function get_pdb_file($job_name) {
	return get_result_file($job_name, ".reconstructed.pdb");
}

function get_output_url($job_name) {
	return "results/".$job_name."/";
}

function print_qstat() {
	global $server_user;
	$cmd = 'qstat -u '.$server_user;
	system($cmd);
}

function get_mb_job_name($job_name) {
	return 'MB-'.$job_name;
}

function is_mb_job_running($job_name) {
	global $server_user;
	$u = escapeshellcmd($server_user);
	$j = escapeshellcmd(get_mb_job_name($job_name));
	$cmd = 'qstat -u '.$u.' -xml | grep "<JB_name>'.$j.'</JB_name>"';
	$result = exec($cmd);
	if (strpos($result, $j) !== false) {
		return true;
	} else {
		return false;
	}
}

if(isset($_GET['job_name'])) {
	$job_name = $_GET['job_name'];
}

if(!isset($job_name)) {
	error("No job name specified");
} else {
	if(file_exists(get_pdb_file($job_name))) {
		header("Location: show_modeling_result.php?job_name=".$job_name);
	}		
}

# auto-refresh
header('Refresh: '.$refresh_rate);
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
<?php
	echo 'Job name='.$job_name.' [<a href="'.get_output_url($job_name).'">see raw output</a>] <br/>';
	
	if(is_mb_job_running($job_name)) {
		echo '<p><b>Your calculation is still running.</b></p><p>This may take a few hours. You can bookmark this page to visit again later. <br/>When the calculation is done the output will appear here. <br/>You may need to refresh the page to see the output.</p>';
	} else {
		echo '<p><font color="red"><b>Your job seems to be not running anymore and probably failed. </b></font></p>';
	}
	echo "Job queue:<br/><br/> (if you don't know what this is, you can safely ignore it)<br/>";
	echo '<p><pre>';
	print_qstat();
	echo '</pre></p>';
?>
<div class="buttons">
<a href="results.php"><img src="images/application_view_list.png" alt="" /> View result list</a>
<a href="index.php"><img src="images/new.png" alt="" /> Start new job</a>
</div>
<div style="clear:both;"></div>
<p></p>
</body>
</html>