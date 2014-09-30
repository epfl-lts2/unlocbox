<?php 
function  div1($pathpage)
{
	if (substr($pathpage,0,10) == substr($_SERVER['PHP_SELF'],0,10))
	{
		echo '<div id="title_selected">';
	}
	else
	{
		echo '<div id="title_unselected">';
	}
}
function  div2($pathpage)
{


		echo '</div>';

}

function print_menu_item($pathpage,$name_item)
{
	div1($pathpage);
	echo '<a href="'.$pathpage.'">'.$name_item.'</a> ';
	div2($pathpage);
}
?>

<div id="topbar"> 
<?php 
print_menu_item("/doc/index.php","Startup");
print_menu_item("/doc/solver/index.php","Solver");
print_menu_item("/doc/prox_operators/index.php","Proximal operators");
print_menu_item("/doc/demos/index.php","Demonstration");
?>
</div> 

