package domain.multiThreading;

import gnu.trove.map.TIntObjectMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;

import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.NucleotideCompound;

import domain.GenomicLocation;
import domain.KmerMap;
import domain.KmerRepeatUnitPair;
import domain.RepUnitBiMap;
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

	public int getID() {
		return threadID;
	}

	@Override
	public void run() {
		System.out.println("THREAD "+threadID+" started on "+ fastq1+"\t"+fastq2);
		try {
			BufferedReader  br1 = new BufferedReader(new FileReader(fastq1));
			BufferedReader  br2 = new BufferedReader(new FileReader(fastq2));
			String line1="";
			String line2="";
			String id="";
			int a = 0;
			int b = 0;
			HashMap<KmerRepeatUnitPair, HashSet<GenomicLocation>> shortKmerMap = generator.getImmidiateFlanksKmerRepUnit2SetOfLocations();
			HashMap<KmerRepeatUnitPair, HashSet<GenomicLocation>> longKmerMap = generator.getLongFlanksKmerRepUnit2SetOfLocations();
			TIntObjectMap<ArrayList<GenomicLocation>> repeatUnitLocationMap = generator.getrepUnitLocationsMap();
			KmerMap kmerIndexMap = generator.getKmerMap();
			RepUnitBiMap repeatIndexingMap = generator.getRepeatMap();
			FindRepetition find_repetition = new FindRepetition();

			int[] allMaxScoresDistribution= new int[1001];
			int[] allUsedScoresDistribution= new int[1001];

			long startTime = System.currentTimeMillis();
			long endTime = System.currentTimeMillis();

			StringTokenizer tokenizer;
			RepetitionList repetitionList = new RepetitionList();

			while ((line1 = br1.readLine())!=null) {

				if (b%50000==0){
					endTime = System.currentTimeMillis();
					System.out.println("THREAD "+threadID+" processed "+b+" lines, (last 50000 lines in "+ (endTime-startTime)/1000 +" seconds");
					startTime = endTime;

				}
				// DEBUGGNG purposes
				/*if (b==20000){
					System.out.println("FINISHED "+b+" LINES");
					break;
				}*/

				line2 = br2.readLine();
				a++;
				b++;

				if (a==1){

					tokenizer = new StringTokenizer(line1.substring(1));
					id = tokenizer.nextToken();
				}
				if (a==2){ // we should examine the sequence
					processRead(line1, line2, id, shortKmerMap, longKmerMap, repeatUnitLocationMap,
							kmerIndexMap, repeatIndexingMap, find_repetition,
							allMaxScoresDistribution,
							allUsedScoresDistribution,repetitionList);
					processRead(line2, line1, id, shortKmerMap, longKmerMap, repeatUnitLocationMap,
							kmerIndexMap, repeatIndexingMap, find_repetition,
							allMaxScoresDistribution,
							allUsedScoresDistribution,repetitionList);
				}
				if (a==4){
					a=0;
				}

			}


		//	bw.close();
		//	fw.close();
			System.out.println("THREAD "+threadID + " finished on "+ fastq1+"\t"+fastq2);

		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void processRead(String line1, String line2, String id,
			TIntObjectMap<HashSet<GenomicLocation>> shortKmerMap,
			TIntObjectMap<HashSet<GenomicLocation>> longKmerMap,			
			TIntObjectMap<ArrayList<GenomicLocation>> repMap,
			KmerMap kmerIndexMap, RepUnitBiMap repeatIndexingMap,
			FindRepetition find_repetition, 
			int[] allMaxScoresDistribution, int[] allUsedScoresDistribution, RepetitionList repetition1)
					throws IOException {
	}
	

}
