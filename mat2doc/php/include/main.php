<?php
// 
function printpage($title,$keywords,$seealso,$demos,$content,$doctype)
{
// define general variable
global $path_include;

// include header
include($path_include."header.php");
print_header($title,$keywords);

//Start of the page
echo'<div id="container">
<div id="header">';
	// top of page
    include($path_include."topofpage.php");
    
 echo'</div>';

    
    

echo '<div id="navigation">';
// top menu
include($path_include."topmenu.php");


echo '</div>';

 echo '<div id="content">';
// print menu
include($path_include."docnav.php");


// -- print the page --
echo '<table width=100% ><tr valign=top><td>'; //create a tabular

//	2) main content
echo '<div id="main_content">';
echo $content;
echo '</div>';
	echo'</td><td>'; // change of cell in the tabular
	
// 	3) right menu (This function-menu)
if ($doctype)
{

	echo'<div id="space"> </div>
	<div id="sidebar">';
	
	echo
	'<b>This function:</b><br>';
	$current_file_name = basename($_SERVER['REQUEST_URI'], ".php"); /* supposing filetype .php*/
	if ($doctype==1)
	{	
		echo '<a href="'.$current_file_name.'_code.php">Program code</a><br>';
	}
	

	if ($doctype==2)
	{	
		echo '<a href="'.substr($current_file_name, 0, strlen($current_file_name)-5).'.php">Help text</a><br>';
	}

	
	if (!empty($seealso))
	{
		echo '<hr />
			<div id="see_also_box">	
			<b>See also:</b> <ul>
			';
		
		foreach($seealso as $cle => $element)
		{
			echo '<li><a href='. $element.'>' . $cle.'</a></li>
			';
		}
		echo '</ul></div>';
	}

	if (!empty($demos))
	{
		echo '<hr />
			<div id="demos_box">	
			<b>Demos:</b> <ul>';
		
		foreach($demos as $cle => $element)
		{
			echo '<li><a href='. $element.'>' . $cle.'</a></li>
			';
		}
		echo '</ul></div>';
	}
	echo '</div>';



}
	
// 	1) left menu (content-menu)
	echo '
	<div id="sidebar">
	<b><u>Subtopics:</u></b><br>
	<hr />';
	include("contentsmenu.php");
	
	$iter=0;
	foreach($menu as $cle => $element)
	{	

		if (substr($cle, 0, 7)=="caption")
		{
			if ($iter==0)
			{
				echo '<b>'.$element.'</b><ul>';	
				$iter=1;
			}
			else
			{
				echo '</ul><b>'.$element.'</b><ul>';
			}
			
		}
		else
		{
			echo '<li>'.$element.'</li>';
		}
	}
		echo '</ul>';
	echo '</div>';
echo'</tr></td></table>'; //close the table
echo '</div>';

// print some footer information
echo '<div id="footer">';
include($path_include."bottomofpage.php");
echo '</div>';


//close the page
echo '</div>';

include($path_include."footer.php");
}
?>






