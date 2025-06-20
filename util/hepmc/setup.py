from util.qol_utils.misc import stdout_redirected
from util.qol_utils.progress_bar import printProgressBar
import subprocess as sub
import sys, os, glob, re, pathlib, importlib
import threading
import atexit
from typing import Optional, Union

# Module-level state for ensuring setup only runs once
_setup_lock = threading.Lock()
_setup_completed = False
_hepmc_available = False
_hepmc_instance = None

# NOTE: pybind11 modules might not play nicely, and once cached we'll
#       encounter issues when uncaching and trying to import a duplicate
#       from another path. I.e. for HepMC3, if we import one from LCG and
#       then try to uncache it and import a different HepMC3 install, we
#       may run into issues. Quite annoying, but it means that this uncaching
#       function isn't really going to be useful. Keeping it here for historical
#       reasons, to remember *why* it doesn't really work.
def uncache_hepmc3():
    modules_to_remove = [key for key in sys.modules.keys() if key.startswith('pyHepMC3')]
    for module in modules_to_remove:
        del sys.modules[module]
        # TODO: Seem to have issues with "import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 as hmroot" after uncaching -- but no issues with "import pyHepMC3". -Jan

def prepend_to_pythonpath(dir:str):
    if(dir is not None and dir not in sys.path):
        sys.path = [dir] + sys.path

def ensure_hepmc_setup(hepmc_dir=None, verbose=False, j=4, require_root=True,
                       require_append=True, require_read_with_index=True, force=False):
    """
    Ensure HepMC3 is set up. Thread-safe and runs only once across all imports.

    Args:
        hepmc_dir: Directory for HepMC installation (None uses default)
        verbose: Enable verbose output
        j: Number of parallel make jobs
        require_root: Require ROOT interface
        require_append: Require append functionality
        require_read_with_index: Require index read functionality
        force: Force reinstallation even if available

    Returns:
        tuple: (success: bool, hepmc_instance: HepMCSetup or None)
    """
    global _setup_completed, _hepmc_available, _hepmc_instance

    if _setup_completed and not force:
        return _hepmc_available, _hepmc_instance

    with _setup_lock:
        if _setup_completed and not force:
            return _hepmc_available, _hepmc_instance

        print("HepMCSetup: Checking HepMC3 availability...")

        # Create the setup instance
        setup_instance = _HepMCSetupInternal(hepmc_dir, verbose)

        try:
            setup_instance.PrepHepMC(
                j=j,
                require_root=require_root,
                require_append=require_append,
                require_read_with_index=require_read_with_index,
                force=force
            )
            _hepmc_available = True
            _hepmc_instance = setup_instance
            print("HepMCSetup: HepMC3 setup completed successfully")
        except Exception as e:
            _hepmc_available = False
            _hepmc_instance = None
            print(f"HepMCSetup: HepMC3 setup failed: {e}", file=sys.stderr)

        _setup_completed = True
        return _hepmc_available, _hepmc_instance


class _HepMCSetupInternal:
    """
    Internal implementation of HepMC setup functionality.
    This class should not be instantiated directly; use HepMCSetup instead.
    """

    def __init__(self, hepmc_dir=None, verbose: bool = False):
        self.version = '3.3.1'  # what version we'll try to install, if necessary
        self.hepmc_dir_default = os.path.dirname(os.path.abspath(__file__)) + '/../../external/hepmc'
        self.python_dir = None
        self.SetVerbose(verbose)

        self.SetDirectory(hepmc_dir)
        self.prefix = "HepMCSetup"

    def SetVerbose(self, verbose: bool):
        self.verbose = verbose

    def SetDirectory(self, hepmc_dir: Optional[str] = None):
        self.hepmc_dir = hepmc_dir
        if self.hepmc_dir is None:
            self.hepmc_dir = self.hepmc_dir_default
        self.hepmc_dir = os.path.normpath(self.hepmc_dir)
        os.makedirs(self.hepmc_dir, exist_ok=True)

        self.logfile = '{}/log.stdout'.format(self.hepmc_dir)
        self.errfile = '{}/log.stderr'.format(self.hepmc_dir)

        self.source_dir = '{}/HepMC3-{}'.format(self.hepmc_dir, self.version)
        self.build_dir = '{}/hepmc3-build'.format(self.hepmc_dir)
        self.install_dir = '{}/hepmc3-install'.format(self.hepmc_dir)

        if self.verbose:
            print('CALLING SetPythonDirectory() within SetDirectory().\tself.install_dir = {}'.format(self.install_dir))
        self.SetPythonDirectory()

    def GetDirectory(self) -> str:
        return self.hepmc_dir

    def SetPythonDirectory(self):
        try:
            self.python_dir = os.path.normpath(glob.glob('{}/**/site-packages'.format(self.install_dir), recursive=True)[0])
        except:
            self.python_dir = None
            if self.verbose:
                print('\tSetPythonDirectory() failed.')
        if self.verbose:
            print('\tSet self.python_dir = {}'.format(self.python_dir))

    def GetPythonDirectory(self) -> Union[str, None]:
        """
        Returns the 'site-packages' directory, where the pyHepMC library is located.
        Adding this to the $PYTHONPATH (sys.path) will allow one to import it.
        """
        return self.python_dir

    def DownloadFork(self):
        """
        Instead of downloading the official HepMC3 release, this will pull
        a (public) fork of the project from CERN GitLab.
        """
        if pathlib.Path(self.source_dir).exists():
            return  # if source exists, don't download again

        with open(self.logfile, 'w') as f, open(self.errfile, 'w') as g:
            # Fetch the HepMC3 source
            hepmc_download = 'https://gitlab.cern.ch/jaofferm/HepMC3/-/archive/jaofferm_writer_append/HepMC3-jaofferm_writer_append.tar.gz'
            hepmc_file = hepmc_download.split('/')[-1]
            if self.verbose:
                self._print('Downloading HepMC3 from {}.'.format(hepmc_download))

            # Depending on Linux/macOS, we use wget or curl.
            has_wget = True
            with stdout_redirected():
                try:
                    sub.check_call('which wget'.split(' '))
                except:
                    has_wget = False
                    pass

            if has_wget:
                sub.check_call(['wget', hepmc_download], cwd=self.hepmc_dir, stdout=f, stderr=g)
            else:
                sub.check_call(['curl', hepmc_download, '-o', hepmc_file], cwd=self.hepmc_dir, stdout=f, stderr=g)
            sub.check_call(['tar', 'xzf', hepmc_file], cwd=self.hepmc_dir, stdout=f, stderr=g)
            sub.check_call(['rm', hepmc_file], cwd=self.hepmc_dir, stdout=f, stderr=g)

            # Rename the source dir to match the one of the official release
            source_dir = '/'.join(self.source_dir.split('/')[:-1]) + '/' + hepmc_file.split('.')[0]
            try:
                sub.check_call(['mv', source_dir, self.source_dir])
            except:  # assume it already exists -- edge case
                sub.check_call(['rm', '-r', source_dir])
                pass

    def Make(self, j: int = 4):
        with open(self.logfile, 'w') as f, open(self.errfile, 'w') as g:

            os.makedirs(self.build_dir, exist_ok=True)

            if self.verbose:
                self._print('Configuring HepMC3 with CMake.')
            # CMake configuration
            command = [
                'cmake',
                '-DCMAKE_INSTALL_PREFIX={}'.format(self.install_dir),
                '-DHEPMC3_ENABLE_ROOTIO:BOOL=ON',
                '-DHEPMC3_ENABLE_PROTOBUFIO:BOOL=OFF',
                '-DHEPMC3_ENABLE_TEST:BOOL=OFF',
                '-DHEPMC3_INSTALL_INTERFACES:BOOL=ON',
                '-DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF',
                '-DHEPMC3_BUILD_DOCS:BOOL=OFF',
                '-DHEPMC3_ENABLE_PYTHON:BOOL=ON',
                '-DHEPMC3_PYTHON_VERSIONS={}.{}'.format(sys.version_info.major, sys.version_info.minor),
                '-DHEPMC3_Python_SITEARCH{a}{b}={c}/lib/python{a}.{b}/site-packages'.format(
                    a=sys.version_info.major, b=sys.version_info.minor, c=self.install_dir),
                self.source_dir
            ]
            for entry in command:
                if entry[0] == '-':
                    self._print('\t{}'.format(entry))
            sub.check_call(command, cwd=self.build_dir, stdout=f, stderr=g)

            # Now make
            if self.verbose:
                self._print('Building HepMC3.')
            command = ['make', '-j{}'.format(j)]
            self.run_make_with_progress(command)

            if self.verbose:
                self._print('Installing HepMC3.')
            command = ['make', 'install']
            sub.check_call(command, cwd=self.build_dir, stdout=f, stderr=g)

    def PrepHepMC(self, j: int = 4, require_root: bool = True, require_append: bool = True,
                  require_read_with_index=True, force=False):
        """
        Checks if pyHepMC3 is available (and has the necessary features),
        and downloads and builds it if not.
        """
        status = True

        # First, check if HepMC Python bindings are available
        try:
            assert self._basic_check()
        except:
            if self.verbose:
                self._print('pyHepMC3 not found.')
            status = False

        # Check ROOT interface if required
        if require_root and status:
            try:
                assert self._root_check()
            except:
                if self.verbose:
                    self._print('pyHepMC3 ROOT interface requested, but not found.')
                status = False

        # Check append functionality if required
        if require_append and status:
            try:
                assert self._root_check_append()
            except:
                if self.verbose:
                    self._print('pyHepMC3 WriterRootTree w/ append functionality requested, but not found.')
                status = False

        # Check read with index functionality if required
        if require_read_with_index and status:
            try:
                assert self._root_check_read()
            except:
                if self.verbose:
                    self._print('pyHepMC3 ReaderRootTree w/ index access functionality requested, but not found.')
                status = False

        if status and not force:
            return

        self.SetDirectory()
        self.DownloadFork()
        self.Make(j)
        self.SetPythonDirectory()

    def _print(self, val: str):
        print('{}: {}'.format(self.prefix, val))

    def _vprint(self, val: str):
        if self.verbose:
            print(val)

    def run_make_with_progress(self, command=['make'], prefix='Building HepMC'):
        # Pattern to match progress indicators
        progress_pattern = re.compile(r'\[\s*(\d+)%\]')

        with open(self.logfile, 'w') as f, open(self.errfile, 'w') as g:
            process = sub.Popen(
                command,
                cwd=self.build_dir,
                stdout=sub.PIPE,
                stderr=sub.PIPE,
                universal_newlines=True,
                bufsize=1  # Line buffered
            )

            # Read stdout line by line
            for line in iter(process.stdout.readline, ''):
                f.write(line)  # Write to log file
                f.flush()

                # Check for progress indicator
                match = progress_pattern.search(line)
                if match:
                    progress = int(match.group(1))
                    printProgressBar(progress, 100, prefix=prefix, suffix='Complete')

            # Read any remaining stderr
            stderr_output = process.stderr.read()
            if stderr_output:
                g.write(stderr_output)

            # Wait for process to complete and check return code
            return_code = process.wait()
            if return_code != 0:
                raise sub.CalledProcessError(return_code, command)

    def _basic_check(self):
        """
        This function checks if (py)HepMC3 is installed
        and ready-to-use.
        """

        if(self.verbose):
            self._vprint('\n' + 10 * '=')
            self._vprint('| BASIC CHECK |')
            self._vprint(10 * '=' + '\n')

            self._vprint('About to uncache HepMC3. Here is sys.path, length {}'.format(len(sys.path)))
            for i,entry in enumerate(sys.path):
                self._vprint('\t[{}] {}'.format(i,entry))

        prepend_to_pythonpath(self.python_dir)

        if(self.verbose):
            self._vprint('Uncached and prepended self.python_dir. Here is sys.path, length {}'.format(len(sys.path)))
            for i,entry in enumerate(sys.path):
                self._vprint('\t[{}] {}'.format(i,entry))
            self._vprint('\t\t(self.python_dir = {} )'.format(self.python_dir))

        # DIAGNOSTIC: Check if the directory exists and what's in it
        if self.python_dir and self.verbose:
            self._vprint("\nDIAGNOSTIC: Checking custom HepMC3 directory: {}".format(self.python_dir))
            if os.path.exists(self.python_dir):
                self._vprint("  Directory exists ✓")
                contents = os.listdir(self.python_dir)
                self._vprint("  Contents: {}".format(contents))

                # Check specifically for pyHepMC3
                pyhepmc3_path = os.path.join(self.python_dir, 'pyHepMC3')
                if os.path.exists(pyhepmc3_path):
                    self._vprint("  pyHepMC3 directory found ✓")
                    pyhepmc3_contents = os.listdir(pyhepmc3_path)
                    self._vprint("  pyHepMC3 contents: {}".format(pyhepmc3_contents))

                    # Check for __init__.py
                    init_file = os.path.join(pyhepmc3_path, '__init__.py')
                    if os.path.exists(init_file):
                        self._vprint("  __init__.py found ✓")
                    else:
                        self._vprint("  __init__.py NOT found ✗")
                else:
                    self._vprint("  pyHepMC3 directory NOT found ✗")
            else:
                self._vprint("  Directory does NOT exist ✗")
                return False

        # DIAGNOSTIC: Try to find pyHepMC3 manually
        self._vprint("\nDIAGNOSTIC: Searching for pyHepMC3 in sys.path...")
        for i, path in enumerate(sys.path):
            potential_path = os.path.join(path, 'pyHepMC3')
            if os.path.exists(potential_path):
                self._vprint("  Found pyHepMC3 at sys.path[{}]: {}".format(i,potential_path))

        # DIAGNOSTIC: Use importlib to find the module
        self._vprint("\nDIAGNOSTIC: Using importlib.util.find_spec...")
        spec = importlib.util.find_spec('pyHepMC3')
        if spec:
            self._vprint("  Module found at: {}".format(spec.origin))
            self._vprint("  Submodule search locations: {}".format(spec.submodule_search_locations))
        else:
            self._vprint("  Module NOT found by importlib ✗")

        try:
            self._vprint("\nDIAGNOSTIC: Attempting import...")
            import pyHepMC3
            self._vprint("  Import successful! File: {}".format(pyHepMC3.__file__))

            from pyHepMC3 import HepMC3 as hm
            self._vprint("  HepMC3 import successful!")

            evt = hm.GenEvent()
            self._vprint("  GenEvent creation successful!")
            return True

        except ImportError as e:
            self._vprint("  Import failed: {}".format(e))
            return False
        except Exception as e:
            self._vprint("  Other error: {}".format(e))
            return False

    def _root_check(self):
        """
        This function checks if the (py)HepMC3 ROOT interface
        is available.
        """

        self._vprint('\n' + 10 * '=')
        self._vprint('| ROOT CHECK |')
        self._vprint(10 * '=' + '\n')

        prepend_to_pythonpath(self.python_dir)

        # DIAGNOSTIC: Check the rootIO structure
        if self.python_dir:
            rootio_path = os.path.join(self.python_dir, 'pyHepMC3', 'rootIO')
            self._vprint("DIAGNOSTIC: Checking rootIO at: {}".format(rootio_path))

            if os.path.exists(rootio_path):
                self._vprint("  rootIO directory exists ✓")
                contents = os.listdir(rootio_path)
                self._vprint("  Contents: {}".format(contents))

                # Check for the specific module
                module_so = os.path.join(rootio_path, 'pyHepMC3rootIO.so')
                module_init = os.path.join(rootio_path, '__init__.py')

                self._vprint("  pyHepMC3rootIO.so exists: {}".format(os.path.exists(module_so)))
                self._vprint("  __init__.py exists: {}".format(os.path.exists(module_init)))

                if os.path.exists(module_so):
                    # Check file permissions and size
                    stat = os.stat(module_so)
                    self._vprint("  .so file size: {} bytes".format(stat.st_size))
                    self._vprint("  .so file permissions: {}".format(oct(stat.st_mode)[-3:]))
            else:
                self._vprint("  rootIO directory does NOT exist ✗")
                return False

        # step-by-step imports -- implicitly doing a _basic_check() again...
        try:
            self._vprint("\nDIAGNOSTIC: Step-by-step import...")

            self._vprint("  Step 1: import pyHepMC3")
            import pyHepMC3
            self._vprint("    ✓ Success: {}".format(pyHepMC3.__file__))

            self._vprint("  Step 2: import pyHepMC3.rootIO")
            import pyHepMC3.rootIO
            self._vprint("    ✓ Success: {}".format(pyHepMC3.rootIO.__file__))

            # For some reason, these checks break on a reload?
            # self._vprint("  Step 3: import pyHepMC3.rootIO.pyHepMC3rootIO")
            # import pyHepMC3.rootIO.pyHepMC3rootIO
            # self._vprint(f"    ✓ Success: {pyHepMC3.rootIO.pyHepMC3rootIO.__file__}")

            # self._vprint("  Step 4: import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3")
            # import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 as hmroot
            # self._vprint("    ✓ Success: HepMC3 imported")

            # self._vprint("  Step 5: access ReaderRootTree")
            # reader_class = hmroot.ReaderRootTree
            # self._vprint("    ✓ Success: ReaderRootTree found")

            return True

        except ImportError as e:
            self._vprint(f"    ✗ ImportError: {e}")

            # Additional diagnostics for ImportError
            import sys
            if hasattr(e, 'name'):
                self._vprint(f"    Failed module: {e.name}")
            if hasattr(e, 'path'):
                self._vprint(f"    Search path: {e.path}")

            return False

        except AttributeError as e:
            self._vprint(f"    ✗ AttributeError: {e}")
            self._vprint("    This suggests the module loaded but doesn't have expected attributes")
            return False

        except Exception as e:
            self._vprint(f"    ✗ Other error: {type(e).__name__}: {e}")
            return False

    def _root_check_append(self):
        """
        This check determines whether or not certain functionality
        is available for the ReaderRootTree/WriterRootTree classes.
        These have been added in a custom fork of HepMC3, which may
        eventually propagate to an official HepMC3 release (and then
        to LCG). Until then, we can make sure to fetch this fork.
        """

        self._vprint('\n' + 10 * '=')
        self._vprint('| ROOT CHECK (APPEND) |')
        self._vprint(10 * '=' + '\n')

        prepend_to_pythonpath(self.python_dir)

        # DIAGNOSTIC: Check the rootIO structure
        if self.python_dir:
            rootio_path = os.path.join(self.python_dir, 'pyHepMC3', 'rootIO')
            self._vprint("DIAGNOSTIC: Checking rootIO at: {}".format(rootio_path))

            if os.path.exists(rootio_path):
                self._vprint("  rootIO directory exists ✓")
                contents = os.listdir(rootio_path)
                self._vprint("  Contents: {}".format(contents))

                # Check for the specific module
                module_so = os.path.join(rootio_path, 'pyHepMC3rootIO.so')
                module_init = os.path.join(rootio_path, '__init__.py')

                self._vprint("  pyHepMC3rootIO.so exists: {}".format(os.path.exists(module_so)))
                self._vprint("  __init__.py exists: {}".format(os.path.exists(module_init)))

                if os.path.exists(module_so):
                    # Check file permissions and size
                    stat = os.stat(module_so)
                    self._vprint("  .so file size: {} bytes".format(stat.st_size))
                    self._vprint("  .so file permissions: {}".format(oct(stat.st_mode)[-3:]))
            else:
                self._vprint("  rootIO directory does NOT exist ✗")
                return False

        # step-by-step imports -- implicitly doing a _basic_check() again...
        try:
            self._vprint("\nDIAGNOSTIC: Step-by-step import...")

            self._vprint("  Step 1: import pyHepMC3")
            import pyHepMC3
            self._vprint("    ✓ Success: {}".format(pyHepMC3.__file__))

            from pyHepMC3 import HepMC3 as hm

            self._vprint("  Step 2: import pyHepMC3.rootIO")
            import pyHepMC3.rootIO
            self._vprint("    ✓ Success: {}".format(pyHepMC3.rootIO.__file__))

            self._vprint("  Step 3: import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3")
            import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 as hmroot
            self._vprint("    ✓ Success.")

            test_file = os.path.dirname(os.path.abspath(__file__)) + '/test/test_events.root'

            self._vprint("  Step 4: Attempt to read event from test file into memory: {}".format(test_file))

            reader = hmroot.ReaderRootTree(test_file)
            evt = hm.GenEvent()
            reader.read_event(evt) # Test file should have 2 events. This will break due to the 2nd argument, if we don't have the right ReaderRootTree::read_event() functionality.
            reader.close()
            self._vprint("    ✓ Success.")

            # copy the test file, we'll append to the copy and then delete it at the end of the test.
            test_file_copy = test_file.replace('test_events.root','test_events_copy.root')
            sub.check_call(['cp',test_file,test_file_copy])

            self._vprint("  Step 5: Attempt to append event to copy of test file: {}".format(test_file_copy))

            writer = hmroot.WriterRootTree(test_file_copy,True) # turns on "append" mode if it's available

            writer.write_event(evt)
            writer.close()
            os.unlink(test_file_copy)
            self._vprint("    ✓ Success.")
            return True

        except ImportError as e:

            try: # clean up just in case
                os.unlink(test_file_copy)
            except:
                pass

            self._vprint(f"    ✗ ImportError: {e}")

            # Additional diagnostics for ImportError
            import sys
            if hasattr(e, 'name'):
                self._vprint(f"    Failed module: {e.name}")
            if hasattr(e, 'path'):
                self._vprint(f"    Search path: {e.path}")

            return False

        except AttributeError as e:

            try: # clean up just in case
                os.unlink(test_file_copy)
            except:
                pass

            self._vprint(f"    ✗ AttributeError: {e}")
            self._vprint("    This suggests the module loaded but doesn't have expected attributes")
            return False

        except Exception as e:

            try: # clean up just in case
                os.unlink(test_file_copy)
            except:
                pass

            self._vprint(f"    ✗ Other error: {type(e).__name__}: {e}")
            return False

    def _root_check_read(self):
        """
        This check determines whether or not certain functionality
        is available for the ReaderRootTree class.
        This has been added in a custom fork of HepMC3, which may
        eventually propagate to an official HepMC3 release (and then
        to LCG). Until then, we can make sure to fetch this fork.
        """

        self._vprint('\n' + 10 * '=')
        self._vprint('| ROOT CHECK (READ) |')
        self._vprint(10 * '=' + '\n')

        prepend_to_pythonpath(self.python_dir)

    # DIAGNOSTIC: Check the rootIO structure
        if self.python_dir:
            rootio_path = os.path.join(self.python_dir, 'pyHepMC3', 'rootIO')
            self._vprint("DIAGNOSTIC: Checking rootIO at: {}".format(rootio_path))

            if os.path.exists(rootio_path):
                self._vprint("  rootIO directory exists ✓")
                contents = os.listdir(rootio_path)
                self._vprint("  Contents: {}".format(contents))

                # Check for the specific module
                module_so = os.path.join(rootio_path, 'pyHepMC3rootIO.so')
                module_init = os.path.join(rootio_path, '__init__.py')

                self._vprint("  pyHepMC3rootIO.so exists: {}".format(os.path.exists(module_so)))
                self._vprint("  __init__.py exists: {}".format(os.path.exists(module_init)))

                if os.path.exists(module_so):
                    # Check file permissions and size
                    stat = os.stat(module_so)
                    self._vprint("  .so file size: {} bytes".format(stat.st_size))
                    self._vprint("  .so file permissions: {}".format(oct(stat.st_mode)[-3:]))
            else:
                self._vprint("  rootIO directory does NOT exist ✗")
                return False

        # step-by-step imports -- implicitly doing a _basic_check() again...
        try:
            self._vprint("\nDIAGNOSTIC: Step-by-step import...")

            self._vprint("  Step 1: import pyHepMC3")
            import pyHepMC3
            self._vprint("    ✓ Success: {}".format(pyHepMC3.__file__))

            from pyHepMC3 import HepMC3 as hm

            self._vprint("  Step 2: import pyHepMC3.rootIO")
            import pyHepMC3.rootIO
            self._vprint("    ✓ Success: {}".format(pyHepMC3.rootIO.__file__))

            self._vprint("  Step 3: import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3")
            import pyHepMC3.rootIO.pyHepMC3rootIO.HepMC3 as hmroot
            self._vprint("    ✓ Success.")

            test_file = os.path.dirname(os.path.abspath(__file__)) + '/test/test_events.root'

            self._vprint("  Step 4: Attempt to read event from test file into memory: {}".format(test_file))

            reader = hmroot.ReaderRootTree(test_file)
            evt = hm.GenEvent()
            reader.read_event(evt) # Test file should have 2 events. This will break due to the 2nd argument, if we don't have the right ReaderRootTree::read_event() functionality.
            reader.close()
            self._vprint("    ✓ Success.")
            return True

        except ImportError as e:
            self._vprint(f"    ✗ ImportError: {e}")

            # Additional diagnostics for ImportError
            import sys
            if hasattr(e, 'name'):
                self._vprint(f"    Failed module: {e.name}")
            if hasattr(e, 'path'):
                self._vprint(f"    Search path: {e.path}")

            return False

        except AttributeError as e:
            self._vprint(f"    ✗ AttributeError: {e}")
            self._vprint("    This suggests the module loaded but doesn't have expected attributes")
            return False

        except Exception as e:
            self._vprint(f"    ✗ Other error: {type(e).__name__}: {e}")
            return False

class HepMCSetup:
    """
    Public interface for HepMC3 setup. This class ensures that expensive
    setup operations only happen once across all imports and instances.

    Usage:
        setup = HepMCSetup()
        if setup.is_available:
            # Use HepMC3
            python_dir = setup.GetPythonDirectory()
    """

    def __init__(self, hepmc_dir=None, verbose: bool = False, j: int = 4,
                 require_root: bool = True, require_append: bool = True,
                 require_read_with_index: bool = True, force: bool = False):
        """
        Initialize HepMCSetup. The actual setup will only run once regardless
        of how many instances are created.

        Args:
            hepmc_dir: Directory for HepMC installation (None uses default)
            verbose: Enable verbose output
            j: Number of parallel make jobs
            require_root: Require ROOT interface
            require_append: Require append functionality
            require_read_with_index: Require index read functionality
            force: Force reinstallation even if available
        """
        self._available, self._instance = ensure_hepmc_setup(
            hepmc_dir=hepmc_dir,
            verbose=verbose,
            j=j,
            require_root=require_root,
            require_append=require_append,
            require_read_with_index=require_read_with_index,
            force=force
        )

    @property
    def is_available(self) -> bool:
        """Returns True if HepMC3 is available and properly set up."""
        return self._available

    def GetPythonDirectory(self) -> Union[str, None]:
        """
        Returns the 'site-packages' directory where pyHepMC is located.
        Adding this to sys.path will allow importing pyHepMC3.
        """
        if self._instance:
            return self._instance.GetPythonDirectory()
        return None

    def GetDirectory(self) -> Union[str, None]:
        """Returns the HepMC installation directory."""
        if self._instance:
            return self._instance.GetDirectory()
        return None

    @staticmethod
    def force_refresh(hepmc_dir=None, verbose=False, j=4, require_root=True,
                      require_append=True, require_read_with_index=True):
        """
        Force a refresh of the HepMC setup. This will re-run all checks
        and potentially reinstall HepMC3.

        Returns:
            HepMCSetup: New instance with refreshed setup
        """
        return HepMCSetup(
            hepmc_dir=hepmc_dir,
            verbose=verbose,
            j=j,
            require_root=require_root,
            require_append=require_append,
            require_read_with_index=require_read_with_index,
            force=True
        )

# Optional: Clean up on exit
atexit.register(lambda: print("HepMCSetup: Cleanup complete") if _setup_completed else None)