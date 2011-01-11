<?php
include_once("inc/settings.php");
session_start();

function is_valid_job_name($job_name) {
	return ereg('^[-0-9a-zA-Z_]+$', $job_name);
}

function is_valid_sequence($sequence) {
	return ereg('^[a-zA-Z]+$', $sequence);
}

function job_output_exists($job_name) {
	return file_exists(get_output_dir($job_name));
}

if (array_key_exists('_submit_check',$_POST)) {

	// save form input
	$seq = trim($_POST['sequence']);
	$seq = str_replace("\n", "", $seq);
	$seq = str_replace("\r", "", $seq);	
	$_SESSION['sequence'] = $_POST['sequence'];
	$_SESSION['job_name'] = substr($_POST['job_name'], 0, 50); # maximum job length is 50

	// check that job name is given
	if(!array_key_exists('job_name',$_POST) || $_POST['job_name'] == "") {
		$form_errors['no_job_name'] = 1;		
	} else {
		// check that job name is valid
		if(!is_valid_job_name($_POST['job_name'])) {
			$form_errors["wrong_job_name"] = 1;
		} else {
			// check that job name is not taken already
			if(job_output_exists($_POST['job_name'])) {
				$form_errors['job_name_exists'] = 1;
			}
		}
	}
	
	// check that sequence exists
	if(!array_key_exists('sequence',$_POST) || $_POST['sequence'] == "") {
		$form_errors["no_sequence"] = 1;
	} else {
		// check that sequence is valid
		if(!is_valid_sequence($seq)) {
			$form_errors['wrong_sequence'] = 1;
		} else {
			$seqlen = strlen($seq);
		}
	}
	
	// check that sequence range is valid
	$seq_from = 1;
	$seq_to = $seqlen;
	// if start or end key exists
	if( (array_key_exists('seqfrom',$_POST) && $_POST['seqfrom'] != "") || 
	    (array_key_exists('seqto',$_POST) && $_POST['seqto'] != "") ){
	    // if start or end key is missing
		if(!array_key_exists('seqfrom',$_POST) || $_POST['seqfrom'] == "" ||
		   !array_key_exists('seqto',$_POST) || $_POST['seqto'] == "") {
			$form_errors['wrong_seq_range'] = 1;
		} else {
			$seq_from = (int) $_POST['seqfrom'];
			$seq_to = (int) $_POST['seqto'];
			if($seq_from < 1 || $seq_to > $seqlen || $seq_to < $seq_from) {
				$form_errors['wrong_seq_range'] = 2;
			}
		}
	} 
	
	// no errors -> redirect to template selection
	if(!isset($form_errors)) {
		$_SESSION['sequence'] = strtoupper($seq);
		$_SESSION['seqfrom'] = $seq_from;
		$_SESSION['seqto'] = $seq_to;
		header("Location: run_template_selection.php");
	}
}
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
  <title>Structure prediction server</title>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
  <link rel="stylesheet" type="text/css" href="css/default.css" />
  
    <script type="text/javascript">
    
    // performed if 'clear form' button is pressed
    function clearForm() {
    	document.seqform.job_name.value='';
    	document.seqform.sequence.value='';
    	document.seqform.seqfrom.value='';
    	document.seqform.seqto.value='';
    	setSeqLen();
    }
    	
	// show the current sequence length in the 'rescount' field
	function setSeqLen() {
		var seqLen = document.getElementById('sequence').value.length;
		var span = document.getElementById('rescount');
		span.innerHTML = seqLen;
	}
	
    // select part of text area
    function selectInTextArea(id, start, end){
    
    	var textarea = document.getElementById(id);
    	if(textarea.setSelectionRange){
        	textarea.setSelectionRange(parseInt(start), parseInt(end));
    	}
    	else {
        	var range = textarea.createTextRange();
        	range.collapse(true);
        
        	range.moveStart('character',parseInt(start) );
        	range.moveEnd('character',(parseInt(end)-parseInt(start)+1));
        	range.select();
    	}
	}
	
	// highlight the subsequence chosen in the sequence range text fields
	function selectSubSequence() {
		var areaId='sequence';
		var fromId = 'seqfrom';
		var toId = 'seqto';
		
		var start = document.getElementById(fromId).value;
		var end = document.getElementById(toId).value;
		
		selectInTextArea(areaId, start-1, end);
	}
	
	// show the start and end indices of the selected sub sequence to the sequence range text fields
	function setRange() {
		var textarea=document.getElementById('sequence');
		document.getElementById('seqfrom').value=textarea.selectionStart+1;
		document.getElementById('seqto').value=textarea.selectionEnd;		
	}
	
	</script>
    
    
</head>
<body onload="setSeqLen();">
<?php include_once("inc/owl_header.inc.php") ?>
<h3>Protein Structure Prediction Server</h3>
<div style="width:800px">
<div class="roundedcornr_box_759934">
   <div class="roundedcornr_top_759934"><div></div></div>
      <div class="roundedcornr_content_759934">

<form name="seqform" action="<?php $_SERVER['PHP_SELF'] ?>" method="POST" enctype="multipart/form-data">
<input type="hidden" name="_submit_check" value="1"/>   
<table cellspacing=10>
  <tbody>
      <tr>
      <td valign="top">Job name:</td>
      <td><input type="text" name="job_name" size=30 maxlength="50" value="<? echo $_SESSION['job_name'] ?>">
      <?php 
      	if(isset($form_errors['no_job_name'])) {
        	echo '<br/><font color="red">Please provide a job name</font>';
      	} 
      	if(isset($form_errors['wrong_job_name'])) {
      		echo '<br/><font color="red">Wrong job name: please use only letters, numbers and underscores</font>';
      	}
      	if(isset($form_errors['job_name_exists'])) {
      		echo '<br/><font color="red">This job name exists already. Please choose a different one.<font>';
      	}      	
      ?>
      </td>
    </tr>
    <tr>
      <td valign="top">Protein sequence:</td> 
      <td><textarea style="width:600px; height:150px" id="sequence" name="sequence" wrap="soft" onselect="setRange();" onkeyup="setSeqLen();" onchange="setSeqLen();"><? echo $_SESSION['sequence'] ?></textarea>
      <?php 
      	if(isset($form_errors['no_sequence'])) {
        	echo '<br/><font color="red" size="2">Please provide a sequence</font>';
      	} 
      	if(isset($form_errors['wrong_sequence'])) {
      		echo '<br/><font color="red" size="2">Illegal character in sequence: please provide a plain protein sequence without fasta header</font>';
      	}      	
      ?>       
      </td>
    </tr>
    <tr>
      <td>Sequence range:</td> 
      <td>
      	<input type="text" name="seqfrom" id="seqfrom" value="<? echo $_SESSION['seqfrom'] ?>" size="4" maxlength="4" onkeyup="selectSubSequence()" onchange="selectSubSequence()" />
      	to <input type="text" name="seqto" id="seqto" value="<? echo $_SESSION['seqto'] ?>" size="4" maxlength="4" onkeyup="selectSubSequence()" onchange="selectSubSequence()" />
      	out of <span id="rescount">0</span>
      	<?php
      	if(isset($form_errors['wrong_seq_range'])) {
        	echo '<br/><font color="red" size="2">Invalid sequence range</font>';
      	} 
      	?>
      </td>
    
    </tr>
    <!-- tr>
      <td>Or upload file:</td>
      <td><INPUT type="file" name="seq_file"></td>
    </tr -->
    <tr>
      <td></td>
      <td><input type="submit" value="Submit" class="fancybutton"> <input type="button" class="fancybutton" value="Clear form" onClick="clearForm();"></td>
    </tr>
  </tbody>
</table>
</form>
</div>
      </div>
   <div class="roundedcornr_bottom_759934"><div></div></div>
</div>
<br/>
<div class="buttons">
<a href="results.php"><img src="images/application_view_list.png" alt="" /> View result list</a>
</div>
<div style="clear:both;"></div>
<p><b>Contact: </b>stehr@molgen.mpg.de</p>
</body>
</html>
