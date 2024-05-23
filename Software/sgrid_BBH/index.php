<?php
 if ($handle = opendir('.')) {
   while (false !== ($file = readdir($handle)))
      {
          //if ($file != "." && $file != "..")
          if ($file != "." && $file != ".." &&
              $file != "index.php" && $file != "index.php~")
          {
                $thelist .= '<a href="'.$file.'">'.$file.'</a>'.'<br>';
          }
       }
  closedir($handle);
  }
?>

<P><a href="..">Go Back</a></p>  

<P>List of files and directories:</p>  
<P><?php echo $thelist;?></p>
