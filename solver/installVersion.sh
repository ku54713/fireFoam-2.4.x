#!/bin/bash
# adds appropriate file and rule to wmake rules directory

destination=$WM_DIR/rules/General/
file=version2

if [ -e  $destination/$file ]
    then
    echo "file $file exists in $destination"
else
    if [ -e $file -a -d $destination ]
	then
	echo "installing $file in $destination"
	cp -rip $file $destination
    fi
fi

file=standard
if [ -e  $destination/$file ]
then
    if ! grep -q "fireFoam" "$destination/$file" ;then
	echo "adding rule to $destination/$file"
	echo "#fireFoam version stamping" >> $destination/$file
	echo "include \$(GENERAL_RULES)/version2" >> $destination/$file
    fi    
fi
