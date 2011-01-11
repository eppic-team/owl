<?php
include_once("inc/settings.php");
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

<font color="red">
<?php
if(isset($_GET['msg'])) {
	$msg = htmlentities(urldecode($_GET['msg']));
	echo 'Internal error: '.$msg;
} else {
	echo 'An unknown internal error occured.';
}
?>
</font>
<p>
<div class="buttons">
<a href="index.php"><img src="images/new.png" alt="" /> Start new job</a>
</div>
<div style="clear:both;"></div>
</div>
</p>
</body>
</html>