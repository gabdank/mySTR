package domain.multiThreading;

import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;

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

import utils.SequenceUtils;
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
				// DEBUGGNG break option
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
					break;
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

	private void processRead(String line1, String mateRead, String id,
			HashMap<KmerRepeatUnitPair, HashSet<GenomicLocation>> shortKmerMap,
			HashMap<KmerRepeatUnitPair, HashSet<GenomicLocation>> longKmerMap,			
			TIntObjectMap<ArrayList<GenomicLocation>> repMap,
			KmerMap kmerIndexMap, RepUnitBiMap repeatIndexingMap,
			FindRepetition find_repetition, 
			int[] allMaxScoresDistribution, int[] allUsedScoresDistribution, RepetitionList repetition1)
					throws IOException {

		System.out.println("processing read1:"+line1);
		repetition1.emptyRepetitionList();
		find_repetition.findMaximalRepetition(line1,repetition1);
		find_repetition.extendMaximalRepetition(id, line1, repetition1);

		// Ensure we are dealing with short tandem repeat (2..8bp long) and that the repeat length is at least 16bp long 
		if (repetition1.get_max_repetition_unit_size()>1 && repetition1.get_max_repetition_unit_size()<9 && repetition1.get_max_repetition_length()>16){
			String repeat1 = line1.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_left_index()+repetition1.get_max_repetition_unit_size());
			if (!repeat1.contains("N")){
				String representativeRepeat1 = SequenceUtils.getRepresentative(repeat1);
				int repeatIndex1 = repeatIndexingMap.getIntForString(representativeRepeat1);
				ArrayList<GenomicLocation> genomicLocations1 = repMap.get(repeatIndex1); 
				if (genomicLocations1 != null){

					String leftFlank = line1.substring(0,repetition1.get_max_repetition_left_index());
					String rightFlank = line1.substring(repetition1.get_max_repetition_right_index());

					String repeatSection = line1.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_right_index() );
					
			
					if (leftFlank.length()>30){
						leftFlank = leftFlank.substring(leftFlank.length()-30, leftFlank.length());
					}
					
					if (rightFlank.length()>30){
						rightFlank = rightFlank.substring(0,30);
					}
					
					/*System.out.println("LEFT FLANK = "+leftFlank);
					System.out.println("RIGHT FLANK = "+rightFlank);
					System.out.println("REPEAT = "+repeatSection);
					*/
					TIntIntMap readFlankKmers = new TIntIntHashMap(); // map that will store all kmers coming from short flanks
					TIntIntMap readPairKmers = new TIntIntHashMap(); // map that will store all kmers from the paired End read

					HashSet<GenomicLocation> kmerUnitPairLocations;

					for (int p=0;p<leftFlank.length()-11;p=p+1){
						String twelveMer = leftFlank.substring(p,p+12);
						if (!twelveMer.contains("N")){
							int kmerIndex = kmerIndexMap.getIntForString(twelveMer);

							kmerUnitPairLocations = shortKmerMap.get(new KmerRepeatUnitPair(kmerIndex, repeatIndex1));

							if (kmerUnitPairLocations != null){										
								if (!readFlankKmers.containsKey(kmerIndex)){
									for (GenomicLocation gl : kmerUnitPairLocations){
										gl.addOneToCounter(threadID);												
									}
									readFlankKmers.put(kmerIndex, 1);
								}
							}

						}
					}

					for (int p=0;p<rightFlank.length()-11;p=p+1){
						String twelveMer = rightFlank.substring(p,p+12);
						if (!twelveMer.contains("N")){	
							int kmerIndex = kmerIndexMap.getIntForString(twelveMer);
							
							kmerUnitPairLocations = shortKmerMap.get(new KmerRepeatUnitPair(kmerIndex, repeatIndex1));
							if (kmerUnitPairLocations != null){							
								if (!readFlankKmers.containsKey(kmerIndex)){
									for (GenomicLocation gl : kmerUnitPairLocations){
										gl.addOneToCounter(threadID);												
									}
									readFlankKmers.put(kmerIndex, 1);
								}
							}
						}
					}


					//--------------------------------------------------------// WORKING HERE!!!					
					// MAXIMUM FINDING				
					double maxScore = -1.0;
					double maxNumOfFlankKmers = readFlankKmers.keySet().size();								
					double maxNumOfPairKmers = readPairKmers.keySet().size();	
					
				
					int numOfLocations = genomicLocations1.size();
					int[] scoresOfLocations=new int[numOfLocations];

					if (maxNumOfFlankKmers>10){// we can rely solely on the short flank kmers - because there are enough of these
						
					}
					else{ // we will have to use the paired end, because the flanks are too short 
						
					}
					
					//--------------------------------------------------------// WORKING HERE
					
					// COUNTER NULLIFICATION
					int[] kmerKeys = readFlankKmers.keys();
					// NULLIFICATION OF THE SHORT FLANK KMERS GENOMIC LOCATIONS
					for (int in = 0;in<kmerKeys.length;in++){									
						int kmerIndex = kmerKeys[in];
						kmerUnitPairLocations = shortKmerMap.get(new KmerRepeatUnitPair(kmerIndex,repeatIndex1));
						if (kmerUnitPairLocations != null){
							for (GenomicLocation gLoc : kmerUnitPairLocations){
								gLoc.nullifyCounter(threadID);
							}
						}
					}
					
					kmerKeys = readPairKmers.keys();
					// NULLIFICATION OF THE PAIRED END READ KMER GENOMIC LOCATIONS
					for (int in = 0;in<kmerKeys.length;in++){									
						int kmerIndex = kmerKeys[in];
						kmerUnitPairLocations = longKmerMap.get(new KmerRepeatUnitPair(kmerIndex,repeatIndex1));
						if (kmerUnitPairLocations != null){
							for (GenomicLocation gLoc : kmerUnitPairLocations){
								gLoc.nullifyCounter(threadID);
							}
						}
					}

				}

			}

		}
	}
}
