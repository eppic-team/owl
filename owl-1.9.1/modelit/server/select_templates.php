<?php
include_once("inc/settings.php");

function get_output_url($job_name) {
	return "results/".$job_name."/";
}

function get_blast_img($job_name) {
	return get_output_url($job_name).$job_name.".blast.png";	
}

function get_psiblast_img($job_name) {
	return get_output_url($job_name).$job_name.".psiblast.png";	
}

function get_blast_img_all($job_name) {
	return get_output_url($job_name).$job_name.".blast.all.png";	
}

function get_psiblast_img_all($job_name) {
	return get_output_url($job_name).$job_name.".psiblast.all.png";	
}

function get_report_file($job_name) {
	return get_result_file($job_name,".report");	
}

function get_blast_file($job_name) {
	return get_result_file($job_name,".pdb-blast.classic.out");	
}

function get_psiblast_file($job_name) {
	return get_result_file($job_name,".pdb-psiblast.classic.out");	
}

function get_template_file($job_name) {
	return get_result_file($job_name, ".man.templates");
}

// takes an array of template ids and writes a template file for modelling
function write_template_file($job_name, $templates) {
	$t_file = get_template_file($job_name);
	$fh = fopen($t_file, 'w') or error('Unable to write template file '.$t_file);	
	foreach($templates as $num => $pdb) {
		fwrite($fh, $pdb."\n");
	}
	fclose($fh);
}

function read_templates($job_name) {
	$result = array();
	if(file_exists(get_template_file($job_name))) {
		$result = file(get_template_file($job_name), FILE_SKIP_EMPTY_LINES | FILE_IGNORE_NEW_LINES);
	}
	return $result;
}

// creates the graphical output from the blast and psiblast output files
function render_blast($job_name) {
	global $num_blast_hits_collapsed;
	global $num_blast_hits_expanded;
	global $eval_cutoff_template_selection;
	$j = escapeshellcmd($job_name);
	$n = $num_blast_hits_collapsed;
	$m = $num_blast_hits_expanded;
	$e = $eval_cutoff_template_selection;
	$w = 900;
	$sc_dir = $GLOBALS['script_dir'];
	$cmd=$sc_dir.'render_blast.pl '.get_blast_file($job_name).' '.$e.' '.$n.' '.$w.' > '.get_blast_img($job_name);
	exec($cmd);
	$cmd=$sc_dir.'render_blast.pl '.get_psiblast_file($job_name).' '.$e.' '.$n.' '.$w.' > '.get_psiblast_img($job_name);
	exec($cmd);
	$cmd=$sc_dir.'render_blast.pl '.get_blast_file($job_name).' '.$e.' '.$m.' '.$w.' > '.get_blast_img_all($job_name);
	exec($cmd);
	$cmd=$sc_dir.'render_blast.pl '.get_psiblast_file($job_name).' '.$e.' '.$m.' '.$w.' > '.get_psiblast_img_all($job_name);
	exec($cmd);
}

if(isset($_GET['job_name'])) {
	$job_name = $_GET['job_name'];
}

if(!isset($job_name)) {
	error("No job name specified");
}

if(isset($_GET['_submit_check'])) {
	// parse templates
	foreach($_GET as $key => $val) {
		if(substr($key,0,4) == "temp" && !empty($val)) {
			$templates[$key] = $val;
		}
	}
	$custom_templates = trim($_GET['custom_templates']);
	if(isset($custom_templates) && !strlen($custom_templates)==0) {
		$cols = split(",",$custom_templates);
		foreach($cols as $col) {
			$temp = trim($col);
			if(!preg_match("/^[0-9][0-9a-z]{3}[A-Z]$/",$temp)) {
				$form_errors['invalid_custom_template'] = strip_tags($temp);			
			} else {
				$templates[] = $temp; # no need to escape shell args because of regexp match above
			}			
		}
	}
	if(!isset($templates) || count($templates) == 0) {
		$form_errors['no_templates'] = 1;
	}
	
	// redirect to run_modeling
	if(!isset($form_errors) && isset($templates) && count($templates) > 0) {
		// write to template file
		write_template_file($job_name, $templates);
		// redirect
		header("Location: run_modeling.php?job_name=".$job_name);
	}
}

render_blast($job_name);

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
  <title>Structure prediction server</title>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
  <link rel="stylesheet" type="text/css" href="css/default.css" />
  
  <script type="text/javascript">
    // preload images
    
    blast_img = new Image();
    blast_img.src = "<? echo get_blast_img($job_name); ?>";

    blast_img_all = new Image();
    blast_img_all.src = "<? echo get_blast_img_all($job_name); ?>";

    psiblast_img = new Image();
    psiblast_img.src = "<? echo get_psiblast_img($job_name); ?>";

    psiblast_img_all = new Image();
    psiblast_img_all.src = "<? echo get_psiblast_img_all($job_name); ?>";
  
    // set up image toggling
  
    var blast_status = true;
    var psiblast_status = true;
    
    function toggle_blast() {
    	blast_status = !blast_status;
    	if(blast_status) {
    		document.getElementById("btn_blast").src = 'images/control_blue_up.png';
    		document.getElementById('img_blast').src = '<? echo get_blast_img_all($job_name); ?>';
    	} else {
    		document.getElementById("btn_blast").src= 'images/control_blue_down.png';
    		document.getElementById('img_blast').src = '<? echo get_blast_img($job_name); ?>';
    	}
    }
    
    function toggle_psiblast() {
    	psiblast_status = !psiblast_status;
    	if(!psiblast_status) {
    		document.getElementById("btn_psiblast").src='images/control_blue_up.png'
    		document.getElementById('img_psiblast').src = '<? echo get_psiblast_img_all($job_name); ?>';
    	} else {
    		document.getElementById("btn_psiblast").src='images/control_blue_down.png'
    		document.getElementById('img_psiblast').src = '<? echo get_psiblast_img($job_name); ?>';
    	}
    }
    
    // show/hide parts of template table
    
    var show_min = <? echo $num_rows_visible_in_template_table; ?>;
    var table_status = true;
    
    function toggle_table() {
    	
    	tbl = document.getElementById("template_table");
		var len = tbl.rows.length;
		    	
    	if(table_status) {
    		document.getElementById("btn_table").src='images/control_blue_up.png'
    		for(i=1 ; i< len; i++){
		 		tbl.rows[i].style.display = "";
	 		} 		
    	} else {
    		document.getElementById("btn_table").src='images/control_blue_down.png'
    		for(i=show_min+1 ; i< len; i++){
		 		tbl.rows[i].style.display = "none";
	 		} 		    		
    	}
    	table_status = !table_status;
    }
    
  </script>
  
</head>
<body>
<?php include_once("inc/owl_header.inc.php") ?>
<div class="main-content">
<?php
	echo 'Job name='.$job_name.' [<a href="'.get_output_url($job_name).'">see raw output</a>] ';
	echo '[<a href="'.get_output_url($job_name).$job_name.'.out.log">output log</a>]<br/>';
?>

<h3>Blast result</h3>
<img id="img_blast" src="<? echo get_blast_img($job_name); ?>" border="0" align="left" />
<a href="Javascript:toggle_blast();"><img id="btn_blast" title="show/hide more hits" src="images/control_blue_down.png" border="0" style="margin-top:10px;" /></a>
<div style="clear:both;"></div>
 
<h3>Psiblast result</h3>
<img id="img_psiblast" src="<? echo get_psiblast_img($job_name); ?>" border="0" align="left" />
<a href="Javascript:toggle_psiblast()"><img id="btn_psiblast" title="show/hide more hits" src="images/control_blue_down.png" border="0" style="margin-top:10px;" /></a>
<div style="clear:both;"></div>
 
<h3>Select templates</h3>
<form action="<?php $_SERVER['PHP_SELF'] ?>" method="GET" enctype="multipart/form-data">
<input type="hidden" name="job_name" value="<? echo $job_name; ?>">
<input type="hidden" name="_submit_check" value="1">
<table border=1 id="template_table" align="left" class="template_table">
  <tbody>

<?php
$lines = file(get_report_file($job_name));

if(file_exists(get_template_file($job_name))) {
	# modeling had already been submitted, show templates being used
	$temps = read_templates($job_name);
	
	if(sizeof($lines) <= 2) {
		echo 'No templates found!<br/><br/>';
//		echo "<p>Templates used for modeling: ";
//		if(count($temps) == 0) {
//			echo "?";
//		} else {
//			foreach($temps as $temp) {
//				echo $temp.' ';
//			}
//			echo "<br/>";
//		}
//		echo '</p>';
		$no_temps = 1;
	} else {	
		for($l = 1; $l < count($lines); $l++) {
			$line = $lines[$l];
			$fields = explode("\t", $line);
			$pdb = $fields[0];
			echo "<tr";
			if($l-1 > $num_rows_visible_in_template_table) echo ' style="display:none" '; 
			echo ">";
			if($l != 1) {
				echo '<td>';
				if(in_array($pdb,$temps)) {
					$key = array_search($pdb, $temps);
					unset($temps[$key]);
					echo '<img src="images/tick.png" alt="yes" />';
				}
				echo '</td>';
			} else {
				echo '<th></th>';
			}
			for($i = 0; $i < count($fields); $i++){
				if($i < 5 || $i > 6) { # skip GTG fields 
					if($l == 1) {
						echo "<th>";
						$field = $fields[$i];
						if(strpos($field, "scop id") !== false) $field = "scop";
						if(strpos($field, "title") !== false) $field = "description";					
						echo $field;
						echo "</th>";
					} else {
						echo "<td>";
						$field = trim($fields[$i]);
						if($i == 0) $field = '<a href="http://pdbwiki.org/index.php/'.substr($pdb,0,4).'">'.$field.'</a>';
						if($i == 7) $field = '<a href="http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sccs='.$field.'">'.$field.'</a>';						
						echo $field;	
						echo "</td>";
					}	
				}
			}
			echo "</tr>\n";
		}	
	}	
	# link to wait for results page
} else {

	# do template selection

	if(sizeof($lines) <= 2) {
		echo 'No templates found!<br/><br/>';
		//echo '<tr><td>Enter a template to use: <input type="text" name="temp1"> e.g. 1tdrA</td></tr>';
		$no_temps = 1;
	} else {
		
		for($l = 1; $l < count($lines); $l++) {
			$line = $lines[$l];
			$fields = explode("\t", $line);
			$pdb = $fields[0];
			echo "<tr";
			if($l-1 > $num_rows_visible_in_template_table) echo ' style="display:none" '; 
			echo ">";
			if($l != 1) {
				echo '<td><INPUT type="checkbox" name="temp'.$l.'" value="'.$pdb.'"></td>';
			} else {
				echo '<th></th>';
			}
			for($i = 0; $i < count($fields); $i++){
				if($i < 5 || $i > 6) { # skip GTG fields 
					if($l == 1) {
						echo "<th>";
						$field = $fields[$i];
						if(strpos($field, "scop id") !== false) $field = "scop";
						if(strpos($field, "title") !== false) $field = "description";						
						echo $field;
						echo "</th>";
					} else {
						echo "<td>";
						$field = trim($fields[$i]);
						if($i == 0) $field = '<a href="http://pdbwiki.org/index.php/'.substr($pdb,0,4).'">'.$field.'</a>';
						if($i == 7) $field = '<a href="http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sccs='.$field.'">'.$field.'</a>';						
						echo $field;
						echo "</td>";
					}
				}
			}
			echo "</tr>\n";
		}
	}	
	
}
?>
  </tbody>
</table><? 
	if(!isset($no_temps)) {
		echo '<a href="Javascript:toggle_table();"><img id="btn_table" title="show/hide more templates" src="images/control_blue_down.png" border="0"/></a>';
		echo '<div style="clear:both;"></div>';	
	}
	?>
<?php if(!file_exists(get_template_file($job_name))) {
	# input user templates
	echo '<table class="template_table"><tbody><tr><td>';
	echo 'Additional templates to use: <input type="text" name="custom_templates"> e.g. 1tdrA, 7dfrA<br/>';	
	echo '</td></tr></tbody></table>';
	
	# error messages
	if(isset($form_errors['invalid_custom_template'])) {
  		echo '<span class="error_msg">Invalid template string. Please enter a comma separated list of pdb codes including chain identifier, e.g. 1tdrA </span><br/>';
  	}
  	if(isset($form_errors['no_templates'])) {
  		echo '<span class="error_msg">Please specify at least one template</span><br/>';  	
  	}
	
	# submit button
	echo '<p><input type="submit" value="Start prediction"></p>';	
} else {
	# output user templates
	echo '<table class="template_table"><tbody><tr><td>';
	echo 'Additional templates used: ';
	if(isset($temps) && count($temps) > 0) {
		foreach($temps as $temp) {
			echo "$temp ";
		}
	} else {
		echo 'none';
	}
	echo '<br/>';	
	echo '</td></tr></tbody></table><br/>';
	
	echo 'Modeling job has been previously submitted with the templates marked above.';
	echo '<div class="buttons">';
	echo '<a href="wait_for_modeling_result.php?job_name='.$job_name.'"><img src="images/eye.png" alt="" /> Check results</a>';
}
?>
</form>
<p>
<div class="buttons">
<a href="index.php"><img src="images/new.png" alt="" /> Start new job</a>
</div>
<div style="clear:both;"></div>
</div>
</body>
</html>
