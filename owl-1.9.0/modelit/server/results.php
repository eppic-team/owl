<?php
include_once("inc/settings.php");

function get_output_url($job_name) {
	return "results/".$job_name."/";
}

function get_seq_file($job_name) {
	return get_result_file($job_name, ".fa");
}

function view_ts_result_url($job) {
	return "select_templates.php?job_name=$job";
}

function view_sm_result_url($job) {
	return "show_modeling_result.php?job_name=$job";
}

function get_yes_img_url() {
	return "images/tick.png";
}

function get_no_img_url() {
	return "images/cross.png";	
}

function get_running_img_url() {
	return "images/anidot.gif";	
}

function get_run_stat_img_url() {
	return "images/stadot.gif";	
}

function get_ts_job_name($job_name) {
	return 'TS-'.$job_name;
}

function is_ts_job_running($job_name) {
	global $server_user;
	$u = escapeshellcmd($server_user);
	$j = escapeshellcmd(get_ts_job_name($job_name));
	$cmd = 'qstat -u '.$u.' -xml | grep "<JB_name>'.$j.'</JB_name>"';
	$result = exec($cmd);
	if (strpos($result, $j) !== false) {
		return true;
	} else {
		return false;
	}
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

# this file should exist if modeling was started
function get_template_file($job_name) {
	return get_result_file($job_name, ".man.templates");
}

# this file should exist if template searching was completed
function get_report_file($job_name) {
	return get_result_file($job_name,".report");	
}

# this file should exist if modelling was completed
function get_pdb_file($job_name) {
	return get_result_file($job_name, ".reconstructed.pdb");
}

# return all job names (= subdirs in results dir) in order of creation date (taken from the fasta file)
function get_all_job_names() {
	$results = array();
	$handler = opendir(get_results_dir());
	while ($file = readdir($handler)) {
        if ($file != '.' && $file != '..' && is_dir(get_output_dir($file))) {
        	$seqfile = get_seq_file($file);
        	if(file_exists($seqfile)) {
        		$last_modified = filemtime($seqfile);
        	} else {
        		$last_modified = 0; 
        	}
        	$results[$file] = $last_modified;
        }
    }
    closedir($handler);
    arsort($results);
    return $results;	
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
<div class="main-content">
<h3>Results overview</h3>
<?php 
	echo '<div class="buttons">';
	echo '<a href="index.php"><img src="images/new.png" alt="" /> Start new job</a>' ;
	echo '</div>';
	echo '<div style="clear:both;"></div>';
	echo '<br/>';

	$jobs = get_all_job_names();
	if(count($jobs) < 1) {
		echo 'No results found!';
	} else {
		echo '<table>';
		echo '<tr>';
        echo '<th>Jobname</th>';
        echo '<th>Submission date</th>';
        echo '<th>Raw output</th>';
        echo '<th>Template selection</th>';
        echo '<th>Modeling result</th>';                
        echo '</tr>';
		foreach($jobs as $job => $last_modified) {
			echo '<tr>';
            
            # job name		
			#echo '<a href=".get_results_url($job).">'.$job.'</a>';
			echo '<td>'.$job.'</td>';
			
			# date
			#$last_modified = filemtime(get_output_dir($job));
			if($last_modified == 0) {
				echo '<td>?</td>';	
			} else {
				echo '<td>'.date("l, dS F, Y @ h:ia", $last_modified).'</td>';				
			}

			
			# raw output
			echo '<td style="text-align:center;"><a href="'.get_output_url($job).'"><img src="images/eye.png" border="0" title="view" alt="view" /></a></td>';
			
			# template selection result
			echo '<td style="text-align:center;">';
			if(is_ts_job_running($job)) {
				echo '<img border="0" src="'.get_running_img_url().'">';
			} else {
				if(!file_exists(get_report_file($job))) {
					echo '<img border="0" src="'.get_no_img_url().'">';
				} else {
					echo '<a href="'.view_ts_result_url($job).'"><img border="0" src="'.get_yes_img_url().'"></a>';
				}					
			}
			echo '</td>';
			
			# modelling result
			echo '<td style="text-align:center;">';
			if(!file_exists(get_template_file($job))) {
				echo '-';
			} else {
				if(is_mb_job_running($job)) {
					echo '<img border="0" src="'.get_running_img_url().'">';
				} else {
					if(!file_exists(get_pdb_file($job))) {
						echo '<img border="0" src="'.get_no_img_url().'">';
					} else {
						echo '<a href="'.view_sm_result_url($job).'"><img border="0" src="'.get_yes_img_url().'"></a>';
					}
				}					
			}		
			echo '</td>';
			echo '</tr>';
		}
		echo '</table>';
		
		echo '<h4>Legend:</h4>';
		echo '<img src="'.get_run_stat_img_url().'"> Job still running<br/>';
		echo '<img src="'.get_yes_img_url().'"> Job completed. Click to see results.<br/>';
		echo '<img src="'.get_no_img_url().'"> No results found. Job probably failed.<br/>';
		echo '<br/>';
		echo 'This page automatically refreshes every '.$refresh_rate.' seconds.';
	}
	
?>
</div>
</body>
</html>
