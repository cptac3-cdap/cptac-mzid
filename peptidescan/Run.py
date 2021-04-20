
from Options import Config, Options, escape_spaces
import subprocess, sys, os, os.path, stat, select, threading, Queue
from Record import Result

class CommandLine(object):
    def __init__(self,options,section):
	self.config = Config(section)
	self.options= options

    def execute(self):
        assert(os.path.exists(self.config.exe))
	if not os.access(self.config.exe,os.X_OK|os.R_OK):
	    os.chmod(self.config.exe,stat.S_IXUSR|stat.S_IXOTH|stat.S_IXOTH|stat.S_IXOTH|stat.S_IXOTH|stat.S_IROTH)
	assert(os.access(self.config.exe,os.X_OK|os.R_OK))
        cmd_string = " ".join([escape_spaces(self.config.exe),
	 	               self.config.opt,
			       str(self.options)])
	ON_POSIX = 'posix' in sys.builtin_module_names
	if self.options.has('v'):
	    print >>sys.stderr, "Executing: %s"%(cmd_string,)
	    sys.stderr.flush()
        process = subprocess.Popen(cmd_string,
	 			   stdin=subprocess.PIPE, 
				   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, 
                                   shell=(sys.platform!="win32"),
				   close_fds=ON_POSIX)
        return process.stdin, process.stdout, process.stderr

def dump(hin,hout):
    for l in hin:
	print >>hout, l.rstrip()
	hout.flush()

def enqueue(handle,queue):
    for line in iter(handle.readline, ''):
        # print >>sys.stderr, line
        queue.put(line)
    handle.close()

class PeptideScanCommandLine(CommandLine):
    def __init__(self,options):
	super(PeptideScanCommandLine,self).__init__(options,"PeptideScan")

    def run(self,peps,blocksize=5000,preprocess=True):

	# Automatically run compress_seq...
	if preprocess:
            csopt = Options()
            csopt.set('i',self.options.get('i'))
	    # csopt.set('I','false')
	    if self.options.has('v'):
		csopt.set('v')
	    if self.options.has('T'):
	        csopt.set('3',36)
	        csopt.set('D','true')
	    else:
	        csopt.set('D','false')
            cs = CompressSeqCommandLine(csopt)
	    cs.run()

	# self.options.set('R','10')

	for b in range(0,len(peps),blocksize):
            qout = Queue.Queue()
	    qerr = Queue.Queue()
	    hin,hout,herr = self.execute()
	    tout = threading.Thread(target=enqueue,args=(hout,qout))
	    tout.daemon = True
            tout.start()
	    terr = threading.Thread(target=enqueue,args=(herr,qerr))
	    terr.daemon = True
            terr.start()
	    hin.write('\n'.join(peps[b:b+blocksize])+'\n')
	    hin.close()
	    
	    while True:
		# stderr
		try:
		    line = qerr.get_nowait()
		    print >>sys.stderr, line.rstrip()
		    continue
		except Queue.Empty:
		    pass
		# stdout
		try:
		    line = qout.get_nowait()
		    yield Result.parse(line)
		    continue
		except Queue.Empty:
		    pass
		# Done?
		if not tout.isAlive() and not terr.isAlive():
		    break
            hout.close()
            herr.close()

class CompressSeqCommandLine(CommandLine):
    def __init__(self,options):
	super(CompressSeqCommandLine,self).__init__(options,"CompressSeq")

    def run(self):
	hin,hout,herr = self.execute()
	hin.close()
	terr = threading.Thread(target=dump,args=(herr,sys.stderr))
	terr.start()
	tout = threading.Thread(target=dump,args=(hout,sys.stdout))
	tout.start()
	terr.join()
	tout.join()
