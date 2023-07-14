* This do-file demonstrates how to install, update, and uninstall the rcr program

****************** INSTALLATION **********************************
* Tell stata to look for ado-files at my website
net from  http://www.sfu.ca/~bkrauth/code
* Describe what's available in the RCR package
net describe rcr
* Install the RCR package
net install rcr.pkg
* Install the ancillary files (i.e., testing and documentation)
net get rcr.pkg, replace

****************** MAINTENANCE **********************************
* ADOUPDATE checks my website to see if the program has been updated since the last download
adoupdate rcr
* If we add the UPDATE option, any changes are applied
adoupdate rcr, update

****************** UNINSTALLATION **********************************
* NET UNINSTALL uninstalls the package
* net uninstall rcr
* The ancillary files are not automatically uninstalled, so delete them by hand.



