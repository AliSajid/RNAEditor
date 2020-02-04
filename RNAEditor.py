#!/usr/bin/python

'''
    Created on May 17, 2013
    Main Class to detect RNA-editing in a given FastQ file
    @author: david
'''

from Helper import Helper, Parameters
from MapFastq import MapFastq
from CallEditingSites import CallEditingSites
import multiprocessing, argparse, os
import traceback
#from PyQt4 import QtGui, QtCore
import textwrap
import sys
import gc
import subprocess


class RnaEdit:

    def __init__(self, fastqFiles, params):
        if isinstance(params, Parameters):
            self.params = params
        else:
            Helper.error("Params has to be Instance of Parameters")


        self.fastqFiles=fastqFiles
        
        #hold the running Popen object
        self.runningCommand=False
        self.isTerminated = False
        #check if the input Files are there
        
        #hold basic statistic values of the run
        basicStatDict={}
            
        
        #set directory where the outputFiles should be written to
        if self.params.output=="default":
            if self.fastqFiles[0].endswith("noDup.realigned.recalibrated.bam"):
                self.sampleName=fastqFiles[0][fastqFiles[0].rfind("/")+1:fastqFiles[0].rfind(".noDup.realigned.recalibrated.bam")]
                self.outdir=fastqFiles[0][0:fastqFiles[0].rfind("/")+1]
            else:
                self.sampleName=fastqFiles[0][fastqFiles[0].rfind("/")+1:fastqFiles[0].rfind(".")]
                # outdir = /path/to/output/rnaEditor/samplename/
                self.outdir=fastqFiles[0][0:fastqFiles[0].rfind("/")+1]+"rnaEditor/"+self.sampleName+"/"
            
            #output=/path/to/output/rnaEditor/samplename/samplename
            self.params.output=self.outdir+self.sampleName
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir, mode=0o755)
                os.chmod(self.outdir, 0o755)

            #create folder for html output
            if not os.path.exists(self.outdir+"/html"):
                os.makedirs(self.outdir +"/html", mode=0o755)
                os.chmod(self.outdir, 0o755)
        
        
        self.checkDependencies()
        
        #check if the input Files are there
        self.printParameters()
        
    def run(self):
        try:
            self.startAnalysis()
        except Exception:
            Helper.error("RnaEditor Failed", self.logFile)
        
        """ At this point the RnaEditor has succesfully finished """
        fileDir = os.path.dirname(os.path.realpath(__file__))
        cmd=["python", fileDir+"/createDiagrams.py", "-o", self.params.output]
        a=subprocess.call(cmd)
        
    def startAnalysis(self):
        """
        START MAPPING
        """
        if self.fastqFiles[0].endswith("bam"):
            if self.fastqFiles[0].endswith("noDup.realigned.recalibrated.bam"):
                Helper.info("Bam File given. Skip mapping", self.logFile)
                self.mapFastQ=None
                mapResultFile=self.fastqFiles[0]
            else: 
                Helper.error("Bam File was not mapped with RnaEditor, this is not supported. Please provide the fastq Files to RnaEditor", self.logFile)
        else:
            self.mapFastQ=MapFastq(self)
            mapResultFile=self.mapFastQ.startAnalysis()

        """
        START CALLING EDITING SITES
        """
        self.callEditSites=CallEditingSites(mapResultFile, self)
        result = self.callEditSites.startAnalysis()
        
        
        
        #finished
        self.isTerminated=True
        
        
        
        Helper.status("rnaEditor Finished with %s" % self.params.output, self.logFile)
        Helper.status("Open %s to see the results" % self.params.output +".html", self.logFile)
        self.cleanUp()
    
    def stopSafely(self):
        Helper.info("Analysis was stopped by User", self.logFile)
    
    def stopImmediately(self):
        if hasattr(self, 'callEditSites'):
            self.callEditSites.cleanUp()
        self.isTerminated=True

        Helper.error("Analysis was terminated by User", self.logFile)
         
    def cleanUp(self):
        #print "deleteAssay " + str(self)
        try:
            if self.mapFastQ!=None:
                self.mapFastQ.cleanUp()
            del self.mapFastQ
        except AttributeError:
            Helper.error("could not delete MapFastQ instance", self.logFile)
        try:
            self.callEditSites.cleanUp()
            del self.callEditSites
        except AttributeError:
            Helper.error("could not delete RnaEdit instance", self.logFile)
        
    def checkDependencies(self):
        """checks if all files are there
        if all programs are installed properly and if the output directory is writable"""
        try:
            self.logFile=open(self.params.output + ".log", "w+")
        except IOError:
            Helper.error("Cannot open Log File")

        if isinstance(self.fastqFiles, list):
            self.fastqFiles=self.fastqFiles
        elif isinstance(self.fastqFiles, str):
            self.fastqFiles=[self.fastqFiles]
        else:
            Helper.error("FastQ File has wrong variable type", self.logFile)
        
        for file in self.fastqFiles:
            if not os.path.isfile(file):
                Helper.error("Could not find: %s" % file, self.logFile)
            
        '''
        Checks the existence of the necessary packages and tools
        :param sourceDir: folder which contains all the software
        '''
        Helper.newline(1)
        Helper.info("CHECK DEPENDENCIES", self.logFile)
        
        #check if all tools are there
        if not os.path.isfile(self.params.sourceDir+"bwa"):
            Helper.error("BWA not found in %s" % self.params.sourceDir, self.logFile)
        if not os.path.isfile(self.params.sourceDir+"picard-tools/SortSam.jar"):
            Helper.error("SortSam.jar not found in %s" % self.params.sourceDir +"picard-tools", self.logFile)
        if not os.path.isfile(self.params.sourceDir+"picard-tools/MarkDuplicates.jar"):
            Helper.error("MarkDuplicates.jar not found in %s" % self.params.sourceDir +"picard-tools", self.logFile)
        if not os.path.isfile(self.params.sourceDir+"GATK/GenomeAnalysisTK.jar"):
            Helper.error("GenomeAnalysisTK.jar not found in %s" % self.params.sourceDir +"GATK/", self.logFile)
        if not os.path.isfile(self.params.sourceDir+"blat"):
            Helper.error("blat not found in %s" % self.params.sourceDir, self.logFile)
        if not os.path.isfile(self.params.sourceDir+"samtools"):
            Helper.error("samtools not found in %s" % self.params.sourceDir, self.logFile)
        if not os.system("java -version")==0:
            Helper.error("Java could not be found, Please install java", self.logFile)
        
        
        
        #check if all files are there
        if not os.path.isfile(self.params.refGenome):
            Helper.error("Could not find Reference Genome in %s: " % self.params.refGenome, self.logFile)
        
        # Files for BWA
        if not os.path.isfile(self.params.refGenome+".amb"):
            Helper.warning("Could not find %s.amb" % self.params.refGenome, self.logFile)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile)
        if not os.path.isfile(self.params.refGenome+".ann"):
            Helper.warning("Could not find %s.ann" % self.params.refGenome, self.logFile)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile)
        if not os.path.isfile(self.params.refGenome+".bwt"):
            Helper.warning("Could not find %s.bwt" % self.params.refGenome, self.logFile)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile)
        if not os.path.isfile(self.params.refGenome+".pac"):
            Helper.warning("Could not find %s.pac" % self.params.refGenome, self.logFile)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile)
        if not os.path.isfile(self.params.refGenome+".sa"):
            Helper.warning("Could not find %s.sa" % self.params.refGenome, self.logFile)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome, self.logFile)

        
        #Files for GATK
        
        
        if self.params.refGenome.endswith("fasta"):
            if not os.path.isfile(self.params.refGenome.replace(".fasta", ".dict")):
                Helper.warning("Could not find %s" % self.params.refGenome.replace(".fasta", ".dict"), self.logFile)
                Helper.error("run: 'java -jar %spicard-tools/CreateSequenceDictionary.jar R=%s  O= %s' to create it" % (self.params.sourceDir, self.params.refGenome, self.params.refGenome.replace(".fastq", ".dict")), self.logFile, self.textField)
        elif self.params.refGenome.endswith("fa"):
            if not os.path.isfile(self.params.refGenome.replace(".fa", ".dict")):
                Helper.warning("Could not find %s" % self.params.refGenome.replace(".fa", ".dict"), self.logFile)
                Helper.error("run: 'java -jar %spicard-tools/CreateSequenceDictionary.jar R=%s  O= %s' to create it" % (self.params.sourceDir, self.params.refGenome, self.params.refGenome.replace(".fa", ".dict")), self.logFile, self.textField)
        else:
            Helper.error("RefGenome has wrong suffix. Either '.fa' or '.fasta'")
        if not os.path.isfile(self.params.refGenome+".fai"):
            Helper.warning("Could not find %s.sai" % self.params.refGenome, self.logFile)
            Helper.error("run: 'samtools faidx %s' to create it" % self.params.refGenome, self.logFile)
    
        #SNP databases
        if not os.path.isfile(self.params.dbsnp):
            Helper.error("Could not find dbSNP database %s: " % self.params.dbsnp, self.logFile)
        if not os.path.isfile(self.params.hapmap) and self.params.hapmap != "None":
            Helper.error("Could not find Hapmap database %s: " % self.params.hapmap, self.logFile)
        if not os.path.isfile(self.params.omni) and self.params.omni != "None":
            Helper.error("Could not find Omni database %s: " % self.params.omni, self.logFile)
        if not os.path.isfile(self.params.esp) and self.params.esp != "None":
            Helper.error("Could not find 1000G database %s: " % self.params.esp, self.logFile)
            
        #region Files
        if not os.path.isfile(self.params.aluRegions):
            Helper.error("Could not find %s: " % self.params.aluRegions, self.logFile)
            
        if not os.path.isfile(self.params.gtfFile):
            Helper.error("Could not find %s: " % self.params.gtfFile, self.logFile)

        
        Helper.info("Dependencies satisfied", self.logFile)
        
        #check if display can activate 

    def printParameters(self):

        Helper.info("*** Start RnaEditor with: ***", self.logFile)
        if self.fastqFiles[0].endswith(".bam"):
            Helper.info("\t Bam File: " + self.fastqFiles[0], self.logFile)
        else:
            if self.params.paired:
                Helper.info("\t FastQ-File_1: " + self.fastqFiles[0], self.logFile)
                Helper.info("\t FastQ-File_2: " + self.fastqFiles[1], self.logFile)
            else:
                Helper.info("\t FastQ-File: " + self.fastqFiles[0], self.logFile)
        Helper.info("\t outfilePrefix:" + self.params.output, self.logFile)
        Helper.info("\t refGenome:" + self.params.refGenome, self.logFile)
        Helper.info("\t dbsnp:" + self.params.dbsnp, self.logFile)
        Helper.info("\t sourceDir:" + self.params.sourceDir, self.logFile)
        Helper.info("\t threads:" + self.params.threads, self.logFile)
        Helper.info("\t maxDiff:" + self.params.maxDiff, self.logFile)
        Helper.info("\t seedDiff:" + self.params.seedDiff, self.logFile)
        Helper.info("\t paired:" + str(self.params.paired), self.logFile)
        Helper.info("\t keepTemp:" + str(self.params.keepTemp), self.logFile)
        Helper.info("\t overwrite:" + str(self.params.overwrite), self.logFile)
        Helper.info("", self.logFile)

       
if __name__ == '__main__':
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(
            prog = 'RnaEditor',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\
             RnaEditor: easily detect editing sites from deep sequencing data'
            ----------------------------------------------------------------
                run without arguments to start the user interface.
            '''))
        parser.add_argument('-i', '--input', metavar='Fastq-Files', nargs='+', type=str, help='Input fastq files (maximum two for paired-end sequencing)', required=True)
        parser.add_argument('-c', '--conf', metavar='Configuration File', type=str, help='Configuration File used to read Parameters for RnaEditor', required=True, default='configuration.txt')
        
        args = parser.parse_args()
        
        parameters = Parameters(args.conf) 
        edit=RnaEdit(args.input, parameters, 0)
        
        edit.start()
        edit.wait()


        del edit
else:
    pass    
    
    