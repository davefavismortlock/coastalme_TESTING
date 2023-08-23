CoastalME (Coastal Modelling Environment) simulates the long-term behaviour of a coast. This initial version is a prototype which considers simple soft cliff cross-shore effects only.

By Andres Payo (British Geological Survey) and Dave Favis-Mortlock.

This is the TESTING version of CoastalME. See <a href="https://github.com/coastalme/coastalme" target="_blank">https://github.com/coastalme/coastalme</a> for the most recent release version of the source code.

INSTRUCTIONS

CoastalME builds easily using Linux. If you wish to run CoastalME on Windows, then we currently recommend using the Cygwin pseudo-Linux software to do this.

1. Create a local copy of the github repository, for example by downloading a zipfile, then unpacking it. We suggest unpacking it to something like "/home/YOUR NAME/coast/CoastalME/", this is your CoastalME folder.

2. In a terminal window (i.e. at a command-line prompt) move to the CoastalME folder (i.e. to "/home/YOUR NAME/coast/CoastalME/" or whatever you decided to call it). Check that you can see a folder structure like this:

      Doxygen<br>
      cshore<br>
      in<br>
      out<br>
      scape<br>
      src<br>

   Then move to the the src folder. Check that you can see a folder structure like this:

      cmake/Modules<br>
      inc<br>
      lib<br>

3. If using Linux: copy run_cmake.sh.LINUX to run_cmake.sh OR if using Cygwin under Windows: copy run_cmake.CYGWIN to run_cmake.sh. Then run <b>run_cmake.sh</b>. If you see error messages re. missing software (for example, telling you that CMake cannot be found or is too old, or GDAL cannot be found or is too old) then you need to install or update the software that is causing the problem.

4. Run <b>make install</b>. This will create an executable file called <b>cme</b> in the CoastalME folder.

5. Edit <b>cme.ini</b> to tell CoastalME which input file to read (for example, in/minimal/minimal.dat)

6. Run <b>cme</b>. Output will appear in the out/ folder.

8. Enjoy!

Dave Favis-Mortlock and Andres Payo

PS We are working on getting CoastalME to build on Windows using the Visual Studio compiler... watch this space!




