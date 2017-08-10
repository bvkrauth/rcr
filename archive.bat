@rem -------------------------------------------------------------
@rem   UPDATE.BAT
@rem   
@rem   Purpose: Archives the active RCR program files, and installs 
@rem            the appropriate files to the web page
@rem   Usage:   ARCHIVE versionnumber
@rem            where "versionnumber" is the desired version number
@rem -------------------------------------------------------------

@rem ----------- Make a new directory ---------
mkdir v.%1
@rem ----------- Copy the files there ---------
copy *.* v.%1

@rem ----------- Update the website ---------
copy *.exe "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code"
copy *.ado "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code"
copy *.hlp "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code"
copy rcr_example.do "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code"
copy rcr_example.dta "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code"
copy rcr.pkg "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code"
copy stata.toc "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code"
echo %1 > "C:\Users\Brian\Documents\My Web Sites\bkrauth-sfu\code\version.txt"

