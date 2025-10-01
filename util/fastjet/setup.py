import os, glob,pathlib, threading, queue, re
from collections import deque
import subprocess as sub
import ROOT as rt
from util.qol_utils.misc import stdout_redirected
from util.qol_utils.progress_bar import printProgressBar, printProgressWithOutput

class FastJetSetup:
    """
    Class for setting up the FastJet library.
    NOTE: We currently download it from GitLab,
    instead of treating it as a submodule.
    """
    def __init__(self,fastjet_dir=None, full_setup=False,verbose=True):
        self.fastjet_version = '3.5.1'
        self.SetDirectory(fastjet_dir)
        if(full_setup): self.PrepFastjet(verbose=verbose)

    def SetDirectory(self,fastjet_dir=None):
        self.fastjet_dir = fastjet_dir
        if(self.fastjet_dir is None):
                self.fastjet_dir = os.path.dirname(os.path.abspath(__file__)) + '/../../external/fastjet'
        self.fastjet_dir = os.path.normpath(self.fastjet_dir)
        os.makedirs(self.fastjet_dir,exist_ok=True)

        self.SetPythonDirectory() # attempts to get the Python library subdir, only sets it if found

        self.logfile = '{}/log.stdout'.format(self.fastjet_dir)
        self.errfile = '{}/log.stderr'.format(self.fastjet_dir)

        self.source_dir = '{}/fastjet-{}'.format(self.fastjet_dir,self.fastjet_version)
        self.build_dir = self.source_dir # a bit weird but it works
        self.install_dir = '{}/fastjet-install'.format(self.fastjet_dir)
        return

    def GetDirectory(self):
        return self.fastjet_dir

    def SetPythonDirectory(self):
        try:
            self.python_dir = os.path.normpath(glob.glob('{}/**/site-packages'.format(self.fastjet_dir),recursive=True)[0])
        except:
            self.python_dir = None

    def GetPythonDirectory(self):
        return self.python_dir

    def PrepFastjet(self, j=4, force=False, verbose=True):
        # Check if Fastjet is already built at destination.
        # Specifically, we will look for some Python-related files.
        if(not force):
            files_to_find = [
                '{}/**/site-packages/fastjet.py'.format(self.fastjet_dir),
                '{}/**/site-packages/_fastjet.a'.format(self.fastjet_dir),
                '{}/**/site-packages/_fastjet.so.0'.format(self.fastjet_dir),
                '{}/**/site-packages/**/*_fastjet*.so*'.format(self.fastjet_dir) # added due to how CMake build works -- not necessary for ./configure method from official release
            ]

            files_to_find = [glob.glob(x,recursive=True) for x in files_to_find]
            files_to_find = [len(x) for x in files_to_find]
            found_fastjet = False
            for entry in files_to_find:
                if(entry > 0):
                    found_fastjet = True # Loose condition since not all these files may exist (e.g. no .so on macos). Should be okay.
                    break
            if(found_fastjet):
                if(verbose): print('Found existing Fastjet installation with Python extension @ {}.'.format(self.fastjet_dir))
                # TODO: Do something else here?
                return

        # Fastjet was not found -> download and build.
        self.DownloadFastjet(verbose=verbose)
        self.BuildFastjet(verbose=verbose)

    def DownloadFastjet(self,force=False,verbose=True):
        with open(self.logfile,'a') as f, open(self.errfile,'a') as g:
            if(not force and pathlib.Path(self.source_dir).exists()): # check if fastjet source is already present, if so we do not need to download it again
                print('Found local fastjet source code (but didn\'t find Python module -- maybe this was not built yet?)')
                return
            # Fetch the Fastjet source
            fastjet_download = 'https://gitlab.com/fastjet/fastjet/-/archive/fastjet-{v}/fastjet-fastjet-{v}.tar.gz'.format(v=self.fastjet_version)
            # fastjet_download = 'http://fastjet.fr/repo/fastjet-{}.tar.gz'.format(self.fastjet_version)
            # fastjet_file = fastjet_download.split('/')[-1]
            fastjet_file = 'fastjet-{}.tar.gz'.format(self.fastjet_version)
            if(verbose): print('Downloading fastjet from {}.'.format(fastjet_download))

            # Depending on Linux/macOS, we use wget or curl.
            has_wget = True
            with stdout_redirected():
                try: sub.check_call('which wget'.split(' '))
                except:
                    has_wget = False
                    pass

            if(has_wget):
                sub.check_call(['wget', fastjet_download,'-O',fastjet_file], cwd=self.fastjet_dir, stdout=f, stderr=g)
            else:
                sub.check_call(['curl',fastjet_download,'-o',fastjet_file], cwd=self.fastjet_dir, stdout=f, stderr=g)

            # Now we force extraction to particular directory name.
            # Unfortunately different Linux distros/macOS have different options
            # available for `tar`, so we'll try different things and hopefully
            # one method works!

            source_dir_short = self.source_dir.replace(self.fastjet_dir + '/','')
            try:
                transform = "s/fastjet-{s}/{s}/".format(s=source_dir_short)
                command = ['tar','zxvf',fastjet_file,'--transform',transform]
                sub.check_call(command, cwd=self.fastjet_dir, stdout=f, stderr=g)
            except:
                os.makedirs(self.source_dir,exist_ok=True)
                command = ['tar','zxvf',fastjet_file,'-C',source_dir_short, '--strip-components=1']
                sub.check_call(command, cwd=self.fastjet_dir, stdout=f, stderr=g)

            sub.check_call(['rm', fastjet_file], cwd=self.fastjet_dir, stdout=f, stderr=g)

    def BuildFastjet(self,j=4,verbose=True):
        """
        To be used with Fastjet downloaded from GitLab,
        which ships with a CMakeLists.txt file.
        Note that we disable SISCone, which is a git submodule
        that we don't get when we download the code as an archive.
        """
        command = [
                'cmake',
                '-DCMAKE_INSTALL_PREFIX={}'.format(self.install_dir),
                '-DFASTJET_ENABLE_PYTHON:BOOL=ON',
                '-DFASTJET_ENABLE_PLUGIN_SISCONE:BOOL=OFF',
                self.source_dir
            ]

        # Note that we request the python bindings in the configuration.
        self.run_command_with_progress(command,cwd=self.source_dir,output_length=45,prefix='Configuring Fastjet:')

        # Now make and install.
        command = ['make', '-j{}'.format(j)]
        self.run_command_with_progress(command,cwd=self.build_dir,prefix='Building Fastjet:',output_width=80)

        command = ['make', 'install']
        self.run_command_with_progress(command,cwd=self.build_dir, output_length=118,prefix='Installing Fastjet:',output_width=80)

        self.SetPythonDirectory()
        return

    def BuildFastjetWithConfigure(self,j=4,verbose=True):
        """
        To be used with the official Fastjet release, which ships
        with a configuration script. (A little old-fashioned, but it
        does get the job done!). This has been deprecated, in favor
        of the CMake-style build.
        """
        # Note that we request the python bindings in the configuration.
        command = ['./configure', '--prefix={}'.format(self.install_dir), '--enable-pyext']
        self.run_command_with_progress(command,cwd=self.source_dir,output_length=324,prefix='Configuring Fastjet:')

        # Now make and install. Will skip "make check".
        self.run_command_with_progress(command=['make', '-j{}'.format(j)],cwd=self.build_dir, output_length=947,prefix='Building Fastjet:',output_width=80)

        self.run_command_with_progress(command=['make', 'install'],cwd=self.build_dir, output_length=458,prefix='Installing Fastjet:',output_width=80)

        self.SetPythonDirectory()
        return

    def run_command_with_progress(self, command, prefix, cwd=None,
                                show_output_lines=10, output_width=80, output_length=None):
        """
        Run command with progress bar and scrolling output display

        Args:
            command: Command to run
            cwd: Current working directory for command
            prefix: Progress bar prefix
            show_output_lines: Number of recent output lines to show (0 to disable)
            output_width: Width of output display area
            output_length: Optional expected length of output; to use if output has no progress tracking (e.g. "[1%] ...")
        """
        # Pattern to match progress indicators
        progress_pattern = re.compile(r'\[\s*(\d+)%\]')

        # Buffer to store recent output lines
        recent_lines = deque(maxlen=show_output_lines * 2)  # Store more than we show for better filtering

        # Track current progress
        current_progress = 0

        def read_stream(stream, file_handle, line_queue, stream_name):
            """Helper function to read from stdout/stderr in separate threads"""
            try:
                for line in iter(stream.readline, ''):
                    if line:
                        file_handle.write(line)
                        file_handle.flush()
                        line_queue.put((stream_name, line))
            except Exception as e:
                line_queue.put(('error', f"Error reading {stream_name}: {str(e)}\n"))
            finally:
                stream.close()

        with open(self.logfile, 'a') as f, open(self.errfile, 'a') as g:

            # add an empty space -- nice to separate outputs from different commands
            f.write('\n')
            g.write('\n')

            process = sub.Popen(
                command,
                cwd=cwd,
                stdout=sub.PIPE,
                stderr=sub.PIPE,
                universal_newlines=True,
                bufsize=1  # Line buffered
            )

            # Create queues for inter-thread communication
            line_queue = queue.Queue()

            # Start threads to read stdout and stderr
            stdout_thread = threading.Thread(
                target=read_stream,
                args=(process.stdout, f, line_queue, 'stdout')
            )
            stderr_thread = threading.Thread(
                target=read_stream,
                args=(process.stderr, g, line_queue, 'stderr')
            )

            stdout_thread.daemon = True
            stderr_thread.daemon = True
            stdout_thread.start()
            stderr_thread.start()

            # Initialize display area
            if show_output_lines > 0:
                # Reserve space for progress bar + output lines
                for i in range(3 + show_output_lines):
                    print()
                # Initial display
                printProgressWithOutput(current_progress, 100, recent_lines,
                                    prefix=prefix, suffix='Complete',
                                    max_lines=show_output_lines, width=output_width)
            else:
                # Just show progress bar
                printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

            counter = 0
            try:
                while process.poll() is None or not line_queue.empty():
                    try:
                        # Get line with timeout to check process status
                        stream_name, line = line_queue.get(timeout=0.1)

                        # Add to recent lines buffer
                        recent_lines.append(line)

                        # Update progress -- either checking for progress indicator,
                        # or expecting some fixed length of output.
                        if(output_length is not None):
                            current_progress = (counter / output_length) * 100.

                        # Check for progress indicator in stdout
                        elif stream_name == 'stdout':
                            match = progress_pattern.search(line)
                            if match:
                                current_progress = int(match.group(1))

                        # Update display (both progress and output)
                        if show_output_lines > 0:
                            printProgressWithOutput(current_progress, 100, recent_lines,
                                                prefix=prefix, suffix='Complete',
                                                max_lines=show_output_lines, width=output_width)
                        else:
                            printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

                        counter += 1
                    except queue.Empty:
                        continue  # Continue checking if process is still running

            except KeyboardInterrupt:
                process.terminate()
                raise

            # Wait for threads to complete
            stdout_thread.join(timeout=1.0)
            stderr_thread.join(timeout=1.0)

            # Final progress update
            if show_output_lines > 0:
                printProgressWithOutput(100, 100, recent_lines,
                                    prefix=prefix, suffix='Complete',
                                    max_lines=show_output_lines, width=output_width)
                print()  # Move cursor to next line after final display
            else:
                printProgressBar(100, 100, prefix=prefix, suffix='Complete')
                print()  # Move to next line

            # Wait for process to complete and check return code
            return_code = process.wait()
            if return_code != 0:
                raise sub.CalledProcessError(return_code, command)

class FastJetInterfaceSetup:
    """
    This class is for setting up the C++/ROOT FastjetInterface code.
    This allows for a custom Fastjet Python interface, via PyROOT,
    that might be faster than using the native Fastjet interface
    (for our current usage patterns).
    """
    def __init__(self,directory=None):
        self.SetDirectory(directory)
        self.executable = '/bin/bash'
        self.lib_extensions = ['so','dylib']
        self.status = False

    def SetDirectory(self,dir):
        if(dir is None):
            dir = os.path.dirname(os.path.abspath(__file__)) + '/../root/fastjet_interface'
        self.dir = dir
        self.logfile = '{}/log.stdout'.format(self.dir)
        self.errfile = '{}/log.stderr'.format(self.dir)

    def FullPreparation(self):
        self.status = False
        try: self.Load(quiet=True)
        except:
            self.Build()
            self.Load()
        self.status = True
        return

    def Build(self):
        """
        Builds the FastjetInterface C++/ROOT library.
        This runs a shell script.
        """
        # TODO: Consider some fancy stuff for setting Pythia8 environment variables?
        self.build_script = 'build.sh'

        command = ['./{}'.format(self.build_script)]
        self.run_command_with_progress(command,cwd=self.dir,prefix='Building FastjetInterface:',output_width=80)
        return

    def Load(self,quiet=False):
        # Load our custom ROOT library.
        try:
            a = rt.FastjetInterface
            return
        except:
            pass

        # Note that we *also* need to fetch some include files -- this has something to do with using the ROOT interpreter.
        # Also note that we have multiple library paths -- to allow for Linux/macOS compatibility.
        build_dir = self.dir + '/fastjet_interface/build'
        include_dir = self.dir + '/fastjet_interface/inc'
        custom_lib_paths = [os.path.realpath('{}/lib/libFastjetInterface.{}'.format(build_dir,x)) for x in self.lib_extensions]
        custom_inc_paths = glob.glob(include_dir + '/**/*.h',recursive=True)

        # Check for any of the libraries.
        found_libary = False
        for libpath in custom_lib_paths:
            if(pathlib.Path(libpath).exists()):
                custom_lib_path = libpath
                found_libary = True
                break

        if(not found_libary):
            if(not quiet): print('Error: The FastjetInterface lib has not been built!')
            assert False

        for inc_path in custom_inc_paths:
            command = '#include "{}"'.format(inc_path)
            status = rt.gInterpreter.Declare(command)
            assert status, 'The following header file did not load properly: {}'.format(inc_path)

        status = rt.gSystem.Load(custom_lib_path)
        assert status == 0, 'The FastjetInterface lib did not load properly.'
        return

    def run_command_with_progress(self, command, prefix, cwd=None,
                                show_output_lines=10, output_width=80, output_length=None):
        """
        Run command with progress bar and scrolling output display

        Args:
            command: Command to run
            cwd: Current working directory for command
            prefix: Progress bar prefix
            show_output_lines: Number of recent output lines to show (0 to disable)
            output_width: Width of output display area
            output_length: Optional expected length of output; to use if output has no progress tracking (e.g. "[1%] ...")
        """
        # Pattern to match progress indicators
        progress_pattern = re.compile(r'\[\s*(\d+)%\]')

        # Buffer to store recent output lines
        recent_lines = deque(maxlen=show_output_lines * 2)  # Store more than we show for better filtering

        # Track current progress
        current_progress = 0

        def read_stream(stream, file_handle, line_queue, stream_name):
            """Helper function to read from stdout/stderr in separate threads"""
            try:
                for line in iter(stream.readline, ''):
                    if line:
                        file_handle.write(line)
                        file_handle.flush()
                        line_queue.put((stream_name, line))
            except Exception as e:
                line_queue.put(('error', f"Error reading {stream_name}: {str(e)}\n"))
            finally:
                stream.close()

        with open(self.logfile, 'a') as f, open(self.errfile, 'a') as g:

            # add an empty space -- nice to separate outputs from different commands
            f.write('\n')
            g.write('\n')

            process = sub.Popen(
                command,
                cwd=cwd,
                stdout=sub.PIPE,
                stderr=sub.PIPE,
                universal_newlines=True,
                bufsize=1  # Line buffered
            )

            # Create queues for inter-thread communication
            line_queue = queue.Queue()

            # Start threads to read stdout and stderr
            stdout_thread = threading.Thread(
                target=read_stream,
                args=(process.stdout, f, line_queue, 'stdout')
            )
            stderr_thread = threading.Thread(
                target=read_stream,
                args=(process.stderr, g, line_queue, 'stderr')
            )

            stdout_thread.daemon = True
            stderr_thread.daemon = True
            stdout_thread.start()
            stderr_thread.start()

            # Initialize display area
            if show_output_lines > 0:
                # Reserve space for progress bar + output lines
                for i in range(3 + show_output_lines):
                    print()
                # Initial display
                printProgressWithOutput(current_progress, 100, recent_lines,
                                    prefix=prefix, suffix='Complete',
                                    max_lines=show_output_lines, width=output_width)
            else:
                # Just show progress bar
                printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

            counter = 0
            try:
                while process.poll() is None or not line_queue.empty():
                    try:
                        # Get line with timeout to check process status
                        stream_name, line = line_queue.get(timeout=0.1)

                        # Add to recent lines buffer
                        recent_lines.append(line)

                        # Update progress -- either checking for progress indicator,
                        # or expecting some fixed length of output.
                        if(output_length is not None):
                            current_progress = (counter / output_length) * 100.

                        # Check for progress indicator in stdout
                        elif stream_name == 'stdout':
                            match = progress_pattern.search(line)
                            if match:
                                current_progress = int(match.group(1))

                        # Update display (both progress and output)
                        if show_output_lines > 0:
                            printProgressWithOutput(current_progress, 100, recent_lines,
                                                prefix=prefix, suffix='Complete',
                                                max_lines=show_output_lines, width=output_width)
                        else:
                            printProgressBar(current_progress, 100, prefix=prefix, suffix='Complete')

                        counter += 1
                    except queue.Empty:
                        continue  # Continue checking if process is still running

            except KeyboardInterrupt:
                process.terminate()
                raise

            # Wait for threads to complete
            stdout_thread.join(timeout=1.0)
            stderr_thread.join(timeout=1.0)

            # Final progress update
            if show_output_lines > 0:
                printProgressWithOutput(100, 100, recent_lines,
                                    prefix=prefix, suffix='Complete',
                                    max_lines=show_output_lines, width=output_width)
                print()  # Move cursor to next line after final display
            else:
                printProgressBar(100, 100, prefix=prefix, suffix='Complete')
                print()  # Move to next line

            # Wait for process to complete and check return code
            return_code = process.wait()
            if return_code != 0:
                raise sub.CalledProcessError(return_code, command)