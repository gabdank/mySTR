package domain.multiThreading;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.NucleotideCompound;

import domain.kmerInformator;
import domain.repeatDetectionMorishita.FindRepetition;
import domain.repeatDetectionMorishita.RepetitionList;

public class processingThread implements Runnable{
	private final int threadID;
	private String fastq1;
	private String fastq2;
	private kmerInformator generator;

	// biojava classes for the alignment
	private final static SubstitutionMatrix<NucleotideCompound> matrix= SubstitutionMatrixHelper.getNuc4_4();
	private static SimpleGapPenalty gapP = new SimpleGapPenalty();

	//morishita's repeat finder
	private final FindRepetition find_repetition = new FindRepetition();
	private final RepetitionList repetition= new RepetitionList();
	
	
	public processingThread(int id,String fastq1File,String fastq2File, kmerInformator generator) throws IOException{
		this.threadID = id;
		this.generator = generator;
		fastq1 = fastq1File;
		fastq2 = fastq2File;	
	
	}
	
	
	@Override
	public void run() {
		System.out.println("Hello from the thread");
		
	}

}
