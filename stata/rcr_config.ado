capture program drop rcr_config
program define rcr_config, rclass
    syntax [, forceos(string) forceversion(string) forcereq(string) nopython notworking]
    di _newline "CONFIGURATION INFORMATION FOR RCR PACKAGE"
    di "------------------------------------------"
    di "The latest RCR package is written in Python. Using it requires " ///
       "that Stata is set"
	di "up to use a compatible version of Python that includes all " ///
	   "required modules."
	di "This command will help ensure that your setup is compatible." _newline
    local version = c(stata_version)
    if "`forceversion'" == "" {
        local version = c(stata_version)
    }
    else {
        local version "`forceversion'"
    }
    if "`forceos'" == "" {
        local os = c(os)
    }
    else {
        local os "`forceos'"
    }
    if "`forcereq'" == "" {
        local requirements "numpy pandas scipy"
    }
    else {
        local requirements "`forcereq'"
    }
    local py = 1

    if `py' == 1 {
        di "STEP 1: Does your Stata version support Python?"
        if (`version' >= 16){
            di "  YES, your Stata version (`version') supports Python." _newline
        }
        else {
            di "  NO, your Stata version (`version') does not support Python."
            di "  Python support is available in Stata versions 16 and above."
			di ""
            local py = 0
        }
    }

    if `py' == 1 {
        di "STEP 2: Do you have Python installed?"
		capture python search
        if _rc == 0 & "`python'" == "" {
            di "  YES, you have Python installed in the following location(s):"
			python search
			di ""
		}
        else {
            di "  NO, Stata is not detecting a Python installation on your " ///
			   "computer."
            di "  You can install one from {browse www.anaconda.com/products/distribution}."
			di "  See {browse blog.stata.com/2020/08/18/stata-python-integration-part-1-setting-up-stata-to-use-python/} "
			di "  for additional instructions and advice." _newline
            local py = 0
        }
	}
    if `py' == 1 {
        di "STEP 3: Is your Stata version linked to a working Python " ///
		   "installation?"
		quietly python query
		local python_exec = r(execpath)
        capture python : 1
        if _rc == 0 & "`tworking'" == "" {
            di "  YES, your Stata version is linked to a working Python " ///
			   "installation at"
            di "  (`python_exec')." _newline
        }
        else {
            di "  NO, your Stata version is not linked to a working Python " ///
			   "installation."
            di "  You can use{help python: python set exec} to use any " ///
			   "installation listed above." _newline
            local py = 0
        }
    }

    if `py' == 1 {
        di "STEP 4: Does your Python installation include all required modules?"
        local missing_modules ""
        foreach module in `requirements' {
            capture python which `module'
            if (_rc > 0) {
                local missing_modules "`module' `missing_modules'"
            }
        }
        if "`missing_modules'" == "" {
            di "  YES, your Python installation includes all required modules."
			di ""
        }
        else {
            di "  NO, your Python installation is missing the following " ///
			   "required modules:"
            di "  `missing_modules'" _newline
			di "To fix this, you have several options:" _newline
			di "Option 1: Use the Anaconda distribution of Python. " ///
			   "Anaconda is designed for use"
			di "          in data analysis, and includes all of the " ///
			   "required modules by default." _newline
			di "          You can download Anaconda at {browse www.anaconda.com/products/distribution}"
			di "          if it is not already installed on your computer." _newline
			di "          You can then use{help python: python set exec} " ///
			   "to tell Stata to use Anaconda."
			di "          (this may require you to restart Stata)" _newline
			di "Option 2: Use the Python Package Installer (pip) to " ///
			   "install the required packages"
			di "          for your current Python installation."
			di "          See {browse blog.stata.com/2020/09/01/stata-python-integration-part-3-how-to-install-python-packages}"
			di "          for instructions and advice." _newline
            local py = 0
        }
    }

    if `py' == 1 {
        di "RCR will use the latest (Python) version." _newline
        di "Add the exe(windows-fortran) or exe(unix-fortran) optional " ///
           "argument if you want "
		di "to use the older (Fortran) version."
        local default_version "python"
    }
    else if (`py' == 0) & "`os'" == "Windows" {
        di "RCR will use the older (Fortran) version."
        local default_version "windows-fortran"
    }
    else if (`py' == 0) & "`os'" == "Unix" {
        di "RCR will use the older (Fortran) version."
        local default_version "unix-fortran"
    }
    else if (`py' == 0) & "`os'" == "MacOSX" {
        di "Unfortunately, the Fortran version of RCR is not available " ///
		   "for the MacOSX"
		di "operating system."
        local default_version "none"
    }
    else if (`py' == 0) {
        di "Unfortunately, the Fortran version of RCR is not available " ///
		   "for the "
		di "`os' operating system."
        local default_version "none"
    }
    return local default_version "`default_version'"
end
