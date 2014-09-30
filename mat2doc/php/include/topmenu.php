<?php 
function  div1m($pathpage)
{
	if ($pathpage == $_SERVER['PHP_SELF'])
	{
		echo '<div id="title_selected2">';
	}
	if ($pathpage == substr($_SERVER['PHP_SELF'],0,4))
	{
		echo '<div id="title_selected2">';
	}
}
function  div2m($pathpage)
{
	if ($pathpage == $_SERVER['PHP_SELF'])
	{
		echo '</div>';
	}
	if ($pathpage == substr($_SERVER['PHP_SELF'],0,4))
	{
		echo '</div>';
	}
}

function print_menu_item2($pathpage,$name_item)
{
	
	echo '<li>';
	div1m($pathpage);
	echo '<a href="'.$pathpage.'">'.$name_item.'</a>';
	div2m($pathpage);
	echo '</li>';
}
?>






<ul>
	<?php 

	print_menu_item2("/index.php","Home");
	print_menu_item2("/doc","Documentation");
	print_menu_item2("/notes/index.php","Notes");
	print_menu_item2("/download/index.php","Download");
	print_menu_item2("mailto:unlocbox-help@lists.sourceforge.net","Contact");
	?>
   <li><a href="http://sourceforge.net/projects/unlocbox">Development page</a></li>
   <li><a href="http://unlocx.math.uni-bremen.de/">UNLocX project page</a></li>
</ul>

