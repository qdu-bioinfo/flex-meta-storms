#!/bin/bash
##Users can change the default environment variables configuration file here
if [[ $SHELL = '/bin/zsh' ]];
then
        PATH_File=~/.zshrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.zsh_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi
else
        PATH_File=~/.bashrc
        if [ ! -f "$PATH_File" ]
        then
                PATH_File=~/.bash_profile
                if [ ! -f "$PATH_File" ]
                then
                        touch $PATH_File
                fi
        fi

fi
PM_PATH=`pwd`
Sys_ver=`uname`
###Checking that environment variable of FlexMetaStorms exists###
Check_old_fms=`grep "export FlexMetaStorms"  $PATH_File|awk -F '=' '{print $1}'`
Check_old_path=`grep "FlexMetaStorms/bin"  $PATH_File |sed 's/\(.\).*/\1/' |awk '{if($1!="#"){print "Ture";}}'`
Check_old_Rscript_path=`grep "FlexMetaStorms/Rscript"  $PATH_File |sed 's/\(.\).*/\1/' |awk '{if($1!="#"){print "Ture";}}'`
#Add_Part="####DisabledbyFlexMetaStorms3####"
echo "**Flex Meta-Storms Installation**"
echo "**version 1.0**"

###Build source code for src package###
if [ -f "Makefile" ]
   then
       echo -e "\n**Flex Meta-Storms src package**"
       make
       echo -e "\n**Build Complete**"
else
   echo -e "\n**Flex Meta-Storms bin package**"
fi
###Configure environment variables###

###Configure environment variables###

if [ "$Check_old_fms" != "" ]
   then
      Checking=`grep ^export\ FlexMetaStorms  $PATH_File|awk -F '=' '{print $2}'`
      if [ "$Checking" != "$PM_PATH" ]
         then
         if [ "$Sys_ver" = "Darwin" ]
            then
            sed -i "" "s/^export\ FlexMetaStorms/$Add_Part\ &/g" $PATH_File
            sed -i "" -e "`grep -n "$Add_Part" $PATH_File | cut -d ":" -f 1 | head -1` a\ 
            export\ FlexMetaStorms=$PM_PATH
            " $PATH_File
         else
             sed -i "s/^export\ FlexMetaStorms/$Add_Part\ &/g" $PATH_File
             sed -i "/$Add_Part\ export\ FlexMetaStorms/a export\ FlexMetaStorms=$PM_PATH" $PATH_File
         fi
     fi    
elif [ "$Check_old_fms" = "" ]
    then
      echo "export FlexMetaStorms="${PM_PATH} >> $PATH_File
fi
if [ "$Check_old_path" = "" ]
    then
      echo "export PATH=\$PATH:\$FlexMetaStorms/bin" >> $PATH_File
fi
if [ "$Check_old_path" = "" ]
    then
	  echo "export PATH=\$PATH:\$FlexMetaStorms/Rscript" >> $PATH_File
fi
###Source the environment variable file###
source $PATH_File
echo -e "\n**Environment Variables Configuration Complete**"
###Configurate the R packages###
echo -e ""
Rscript $PM_PATH/Rscript/config.R
###End
echo -e "\n**Flex Meta-Storms Installation Complete**"
echo -e "\n**An example dataset with demo script is available in \"example\"**"
