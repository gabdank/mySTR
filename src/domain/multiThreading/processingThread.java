package domain.multiThreading;

import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.THashSet;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashBigSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;

import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;

import utils.SequenceUtils;
import domain.AlignmentsData;
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
	private FileWriter fw;
	private BufferedWriter bw;

	public processingThread(int id,String fastq1File,String fastq2File, kmerInformator generator, String pathOutput) throws IOException{
		this.threadID = id;
		this.generator = generator;
		fastq1 = fastq1File;
		fastq2 = fastq2File;	
		fw = new FileWriter(new File(pathOutput+id+".ou"));
		bw = new BufferedWriter(fw); 
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
			Object2ObjectOpenHashMap<KmerRepeatUnitPair, ObjectOpenHashBigSet<GenomicLocation>> shortKmerMap = generator.getImmidiateFlanksKmerRepUnit2SetOfLocations();
			Object2ObjectOpenHashMap<KmerRepeatUnitPair, ObjectOpenHashBigSet<GenomicLocation>> longKmerMap = generator.getLongFlanksKmerRepUnit2SetOfLocations();
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
					//System.out.println("Size of ");
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
					//if (id.equals("chr1.AG-3248/1")){
						processRead(line1, line2, id, shortKmerMap, longKmerMap, repeatUnitLocationMap,
							kmerIndexMap, repeatIndexingMap, find_repetition,
							allMaxScoresDistribution,
							allUsedScoresDistribution,repetitionList);
						processRead(line2, line1, id, shortKmerMap, longKmerMap, repeatUnitLocationMap,
							kmerIndexMap, repeatIndexingMap, find_repetition,
							allMaxScoresDistribution,
							allUsedScoresDistribution,repetitionList);
					//}
					//break;
				}
				if (a==4){
					a=0;
				}

			}


			bw.close();
			fw.close();
			System.out.println("THREAD "+threadID + " finished on "+ fastq1+"\t"+fastq2);

		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void processRead(String line1, String mateRead, String id,
			Object2ObjectOpenHashMap<KmerRepeatUnitPair, ObjectOpenHashBigSet<GenomicLocation>> shortKmerMap,
			Object2ObjectOpenHashMap<KmerRepeatUnitPair, ObjectOpenHashBigSet<GenomicLocation>> longKmerMap,			
			TIntObjectMap<ArrayList<GenomicLocation>> repMap,
			KmerMap kmerIndexMap, RepUnitBiMap repeatIndexingMap,
			FindRepetition find_repetition, 
			int[] allMaxScoresDistribution, int[] allUsedScoresDistribution, RepetitionList repetition1)
					throws IOException {

		//System.out.println("processing read1:"+line1);
		repetition1.emptyRepetitionList();
		find_repetition.findMaximalRepetition(line1,repetition1);
		
		find_repetition.extendMaximalRepetition(id, line1, repetition1);
		//System.out.println("unit size : "+repetition1.get_max_repetition_unit_size());
		//System.out.println("length of the repetition : "+repetition1.get_max_repetition_length());
		// Ensure we are dealing with short tandem repeat (2..8bp long) and that the repeat length is at least 16bp long 
		if (repetition1.get_max_repetition_unit_size()>1 && repetition1.get_max_repetition_unit_size()<9 && repetition1.get_max_repetition_length()>16){
			String repeat1 = line1.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_left_index()+repetition1.get_max_repetition_unit_size());
			//System.out.println("Found repetition:"+repeat1);
			if (!repeat1.contains("N")){
				//System.out.println("No Ns in the repeat :)");
				String representativeRepeat1 = SequenceUtils.getRepresentative(repeat1);
				
				int repeatIndex1 = repeatIndexingMap.getIntForString(representativeRepeat1);
				ArrayList<GenomicLocation> genomicLocations1 = repMap.get(repeatIndex1); 
				if (genomicLocations1 != null){
					//System.out.println("There are some potential genomic locatons");
					String leftFlank = line1.substring(0,repetition1.get_max_repetition_left_index());
					String rightFlank = line1.substring(repetition1.get_max_repetition_right_index());

					String repeatSection = line1.substring(repetition1.get_max_repetition_left_index(),repetition1.get_max_repetition_right_index() );


					if (leftFlank.length()>30){
						leftFlank = leftFlank.substring(leftFlank.length()-30, leftFlank.length());
					}

					if (rightFlank.length()>30){
						rightFlank = rightFlank.substring(0,30);
					}

					//System.out.println("LEFT FLANK = "+leftFlank);
					//System.out.println("RIGHT FLANK = "+rightFlank);
					//System.out.println("REPEAT = "+repeatSection);
					
					TIntIntMap readFlankKmers = new TIntIntHashMap(); // map that will store all kmers coming from short flanks
					TIntIntMap readPairKmers = new TIntIntHashMap(); // map that will store all kmers from the paired End read

					ObjectOpenHashBigSet<GenomicLocation> kmerUnitPairLocations;

					// Going across left flank and extracting kmers from it, creating pairs with the repeatUnit and retreiving the genomic locations
					// matching these pairs, adding 1 to each genomic location counter
					// and record the kmers that were used for this read
					
					for (int p=0;p<leftFlank.length()-11;p=p+1){
						String twelveMer = leftFlank.substring(p,p+12);
						//System.out.println("left flank Kmer "+twelveMer);
						if (!twelveMer.contains("N")){
							int kmerIndex = kmerIndexMap.getIntForString(twelveMer);
							//System.out.println("left flank Kmer index "+kmerIndex);
							kmerUnitPairLocations = shortKmerMap.get(new KmerRepeatUnitPair(kmerIndex, repeatIndex1));
							//System.out.println("Pair = " + new KmerRepeatUnitPair(kmerIndex, repeatIndex1));
							
							
							
							if (kmerUnitPairLocations != null){
								//System.out.println("Num of locatoins of the pair "+kmerUnitPairLocations.size());
								if (!readFlankKmers.containsKey(kmerIndex)){
									for (GenomicLocation gl : kmerUnitPairLocations){
										gl.addOneToCounter(threadID);												
									}
									readFlankKmers.put(kmerIndex, 1);
								}
							}

						}
					}
					//System.out.println("number of flank kmers detected was "+readFlankKmers.size());
					// same story with right flank
					
					for (int p=0;p<rightFlank.length()-11;p=p+1){
						String twelveMer = rightFlank.substring(p,p+12);
						//System.out.println("right flank Kmer "+twelveMer);

						if (!twelveMer.contains("N")){	
							int kmerIndex = kmerIndexMap.getIntForString(twelveMer);
							//System.out.println("Pair = " + new KmerRepeatUnitPair(kmerIndex, repeatIndex1));

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
					//System.out.println("number of flank kmers detected was "+readFlankKmers.size());

										
					// MAXIMUM FINDING	
					// Now that we have readFlankMers set, we would like to find potential locations for alignment
					
					double maxScore = -1.0;
					double maxNumOfFlankKmers = readFlankKmers.keySet().size();								
					double maxNumOfPairKmers = readPairKmers.keySet().size();	


					int numOfLocations = genomicLocations1.size();
					int[] scoresOfLocations=new int[numOfLocations];
					
					//System.out.println("Number of k-mers was:"+maxNumOfFlankKmers);
					if (maxNumOfFlankKmers>15){// we can relay solely on the short flank kmers - because there are enough of these
						int ind = 0;
						for (GenomicLocation gl : genomicLocations1){
							int localKmerCount = gl.getCounter(threadID);
							if (maxScore<localKmerCount){
								maxScore = localKmerCount;
							}			
							scoresOfLocations[ind]=localKmerCount;
							ind++;
						}	
						
						// finding the maximum score in the genomic locations
						// adding "close" to the maximum value scores to the list of potential locations
						ArrayList<GenomicLocation> potentialMappings = new ArrayList<GenomicLocation>();
						for (GenomicLocation gl : genomicLocations1){
							int localKmerCount = gl.getCounter(threadID);
							if (localKmerCount>=maxScore-2){
								potentialMappings.add(gl);
							}
						}

						
						detectAndReportAlignment(potentialMappings, id, line1, mateRead,leftFlank, rightFlank, repeatSection);

					}
					else{ // we will have to use the paired end, because the flanks are too short 
						// rewrite the solution for usage of paired end.
						// ideally we should go for sparse count of kmers - let say every 5bp or something like - that 
						// ensuring fast count and selection of potential Mappings
						// we should not forget the kmers counts coming from the flanks - they are not many - but still there
						// we shouldn't attemp aligning without taking into consideration the short flanks sequence kmers
						//System.out.println("Going through the paired end - since flanks were not enough informative");
						
						for (int p=0;p<mateRead.length()-11;p=p+1){
							String twelveMer = mateRead.substring(p,p+12);
							if (!twelveMer.contains("N")){
								int kmerIndex = kmerIndexMap.getIntForString(twelveMer);
								kmerUnitPairLocations = longKmerMap.get(new KmerRepeatUnitPair(kmerIndex, repeatIndex1));
								if (kmerUnitPairLocations!=null){
									if (!readPairKmers.containsKey(kmerIndex)){
										for (GenomicLocation gl : kmerUnitPairLocations){
											gl.addOneToCounter(threadID);												
										}
										readPairKmers.put(kmerIndex, 1);
									}
								}
							}
						}
						maxNumOfPairKmers = readPairKmers.keySet().size();
						
						//System.out.println("Maximum number of kmers for the pair end: "+maxNumOfPairKmers);
						
						maxScore = -1.0;
						if (maxNumOfFlankKmers+maxNumOfPairKmers>20){// we can relay on mixture of flank Kmers and paired end kmers
							int ind = 0;
							for (GenomicLocation gl : genomicLocations1){
								int localKmerCount = gl.getCounter(threadID);
								if (maxScore<localKmerCount){
									maxScore = localKmerCount;
								}			
								scoresOfLocations[ind]=localKmerCount;
								ind++;
							}	
							
							// finding the maximum score in the genomic locations
							// adding "close" to the maximum value scores to the list of potential locations
							ArrayList<GenomicLocation> potentialMappings = new ArrayList<GenomicLocation>();
							for (GenomicLocation gl : genomicLocations1){
								int localKmerCount = gl.getCounter(threadID);
								if (localKmerCount>=maxScore-2){
									potentialMappings.add(gl);
								}
							}

							
							detectAndReportAlignment(potentialMappings, id, line1, mateRead,leftFlank, rightFlank, repeatSection);
						}
					}


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

				} // genomic location was NULL end of the 

			}

		}
		
	}

	/** 
	 * Make sure the length of read flanks if not bigger than the genomic location flanks - otherwise we can not use them as anchors
	 * 
	 */
	private boolean validLengths(String readLeft, String readRight, String glLeft, String glRight){
		int readMAximalFlanksLength = Math.max(readLeft.length(), readRight.length());
		int glFlank = Math.min(glLeft.length(),glRight.length());
		return readMAximalFlanksLength<=glFlank;
		
	}
	/**
	 * The method will take a read pair and ArrayList of genomic locations it will attept to align the reads to
	 * in case promising alignment found it will print the genomic location and the alignment to the file created by the running thread.
	 * These files will be later merged.
	 * 
	 * For the alignment the following data is needed:
	 * 	data[i]=align(readsLeftFlank, readsRightFlank, readsRepetitiveReadSeq, genomicLocation.getLeftFlank(), 
	 * genomicLocation.getRightFlank(), genomicLocation.getRepetitiveSequence(), readsPairedEnd, matrix, gapP,readID);
	 * We can elevate the gap penalty on the extension - to prevent gapped alignment
	 * Other than that I am not sure
	 * @throws IOException 
	 */
	private void detectAndReportAlignment(ArrayList<GenomicLocation> potentialLocations, String readID, String read, 
											String pairRead,String leftFlank, String rightFlank, String repetitiveSection) throws IOException{
		// if there is high match print into FileWriter
		// probably print something like what we have done earlier for Refined STRExpansion
		ArrayList<AlignmentsData> allignmentDataList = new ArrayList<AlignmentsData>();
		
		for (GenomicLocation gl : potentialLocations){
			// TODO align and get the alignment score
			// TODO decide later what to do with it
			if (validLengths(leftFlank, rightFlank,gl.getLeftFlank(), gl.getRightFlank())==true){
				//System.out.println(align(leftFlank, rightFlank, repetitiveSection, gl.getLeftFlank(), gl.getRightFlank(), gl.getRepetitiveSection(), pairRead, matrix,	gapP, readID));
				allignmentDataList.add(align(leftFlank, rightFlank, repetitiveSection, gl.getLeftFlank(), gl.getRightFlank(), gl.getRepetitiveSection(), pairRead, matrix,	gapP, readID));
			}
		}

		// Finding the maximum score and adding the alignment to the list of the corresponding genomic location
		double maximalScore = 0;
		int maximalIndex = -1;

		for (int ind=0;ind<allignmentDataList.size(); ind++){
			if (maximalScore<allignmentDataList.get(ind).getScore()){
				maximalScore = allignmentDataList.get(ind).getScore();
				maximalIndex = ind;
			}								
		}
		// check if there are no 2 almost exact maximums
		int maximumsCounter = 0;
		for (int ind=0;ind<allignmentDataList.size();ind++){
			if (Math.abs(maximalScore-allignmentDataList.get(ind).getScore())<0.00001){
				maximumsCounter++;
			}
		}
		
		//System.out.println("Number of maximums:"+maximumsCounter);
		//System.out.println("Maximal score = "+maximalScore);
		//System.out.println("------------------");
		
		if (maximumsCounter==1 && maximalScore>=0.85){			
			bw.write(potentialLocations.get(maximalIndex).reportTheGenomicLocation()+"{[");
			bw.write(allignmentDataList.get(maximalIndex).reportTheAlignmentData()+"]#}\n");							
		}
		

	}
	
	private static AlignmentsData align(String readLeftFlank, String readRightFlank, String readRepetitiveRegion, String genomicLeftFlank, 
			String genomicRightFlank, String genomicRepetitiveRegion, String pairedEnd, SubstitutionMatrix<NucleotideCompound> mat,	SimpleGapPenalty gap, String id){
		// forward alignment
		int rLeftFlankLen = readLeftFlank.length();
		int rRightFlankLen = readRightFlank.length();
		//System.out.println("STARTING TO ALIGN:");
		//System.out.println();
		//System.out.println(readLeftFlank+" Left length "+rLeftFlankLen);
		//System.out.println(readRightFlank+" Right length "+rRightFlankLen);
					
		DNASequence readLeft = new DNASequence(readLeftFlank+readRepetitiveRegion.substring(0,8), AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence readRight = new DNASequence(readRepetitiveRegion.substring(readRepetitiveRegion.length()-8)+readRightFlank, AmbiguityDNACompoundSet.getDNACompoundSet());
		
		//int genomicLeftLen = readLeft.getLength();
		//int genomicRightLen = readRight.getLength();
		
			
		DNASequence genomicLeft = new DNASequence(genomicLeftFlank.substring(genomicLeftFlank.length()-rLeftFlankLen)+genomicRepetitiveRegion.substring(0,8), AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence genomicRight = new DNASequence(genomicRepetitiveRegion.substring(genomicRepetitiveRegion.length()-8)+genomicRightFlank.substring(0,rRightFlankLen), AmbiguityDNACompoundSet.getDNACompoundSet());
		
		SequencePair<DNASequence, NucleotideCompound> leftSideForward =Alignments.getPairwiseAlignment(readLeft, genomicLeft,PairwiseSequenceAlignerType.GLOBAL, gap,mat);
		SequencePair<DNASequence, NucleotideCompound> rightSideForward =Alignments.getPairwiseAlignment(readRight, genomicRight,PairwiseSequenceAlignerType.GLOBAL, gap,mat);
		
		//System.out.println("left fwd :\n"+leftSideForward);
		//System.out.println("right fwd :\n"+rightSideForward);
		//reverse alignment
		String reverseReadLeftSeq = SequenceUtils.getReverseComplementary(readLeftFlank);
		String reverseReadRepetitiveSeq = SequenceUtils.getReverseComplementary(readRepetitiveRegion);
		String reverseReadRightSeq = SequenceUtils.getReverseComplementary(readRightFlank);

		DNASequence reverseReadLeft = new DNASequence(reverseReadRightSeq+reverseReadRepetitiveSeq.substring(0,8), AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence reverseReadRight = new DNASequence(reverseReadRepetitiveSeq.substring(reverseReadRepetitiveSeq.length()-8)+reverseReadLeftSeq, AmbiguityDNACompoundSet.getDNACompoundSet());

		DNASequence revGenomicLeft = new DNASequence(genomicLeftFlank.substring(genomicLeftFlank.length()-rRightFlankLen)+genomicRepetitiveRegion.substring(0,8), AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence revGenomicRight = new DNASequence(genomicRepetitiveRegion.substring(genomicRepetitiveRegion.length()-8)+genomicRightFlank.substring(0,rLeftFlankLen), AmbiguityDNACompoundSet.getDNACompoundSet());

		SequencePair<DNASequence, NucleotideCompound> leftSideReverse =Alignments.getPairwiseAlignment(reverseReadLeft, revGenomicLeft,PairwiseSequenceAlignerType.GLOBAL, gap,mat);
		SequencePair<DNASequence, NucleotideCompound> rightSideReverse =Alignments.getPairwiseAlignment(reverseReadRight, revGenomicRight,PairwiseSequenceAlignerType.GLOBAL, gap,mat);

		// Align the paired end to the genomic sequence with LOCAL alignment
		DNASequence pairedRead = new DNASequence(pairedEnd, AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence reversePairedRead = new DNASequence(SequenceUtils.getReverseComplementary(pairedEnd), AmbiguityDNACompoundSet.getDNACompoundSet());

		DNASequence genomicLeftFull = new DNASequence(genomicLeftFlank, AmbiguityDNACompoundSet.getDNACompoundSet());
		DNASequence genomicRightFull = new DNASequence(genomicRightFlank, AmbiguityDNACompoundSet.getDNACompoundSet());

		SequencePair<DNASequence, NucleotideCompound> pairedLeft =Alignments.getPairwiseAlignment(pairedRead, genomicLeftFull,PairwiseSequenceAlignerType.LOCAL, gap,mat);
		SequencePair<DNASequence, NucleotideCompound> pairedRight =Alignments.getPairwiseAlignment(reversePairedRead, genomicRightFull,PairwiseSequenceAlignerType.LOCAL, gap,mat);

		//System.out.println("Left pair \n"+pairedLeft);
		//System.out.println((double)pairedLeft.getNumIdenticals()/pairedEnd.length());
		//System.out.println("Right pair \n"+pairedRight);
		//System.out.println((double)pairedRight.getNumIdenticals()/pairedEnd.length());
		/*
		 * assuming that everything worked fine we can now evaluate our alignments
		 * so there are two main cases - forward flank and reverse paired end or reversed flanks and forward paired end
		 */
		
		int totalForwardFlanksLength=-1;
		int totalForwardFlanksIdenticals=-1;
		int totalReverseFlanksLength=-1;
		int totalReverseFlanksIdenticals=-1;
		
		int flankSizeLimit = 5;
		boolean noFlanks = false;
		
		if (rLeftFlankLen>=flankSizeLimit && rRightFlankLen>=flankSizeLimit){
			totalForwardFlanksLength = leftSideForward.getLength()+ rightSideForward.getLength();
			totalForwardFlanksIdenticals = leftSideForward.getNumIdenticals()+ rightSideForward.getNumIdenticals();
			totalReverseFlanksLength = leftSideReverse.getLength()+ rightSideReverse.getLength();
			totalReverseFlanksIdenticals = leftSideReverse.getNumIdenticals()+ rightSideReverse.getNumIdenticals();
		}
		else{
			if (rLeftFlankLen>=flankSizeLimit && rRightFlankLen<flankSizeLimit){
				//System.out.println("left is ok, right is small");
				
				totalForwardFlanksLength = leftSideForward.getLength();
				totalForwardFlanksIdenticals = leftSideForward.getNumIdenticals();
					
				totalReverseFlanksLength = rightSideReverse.getLength();
				totalReverseFlanksIdenticals = rightSideReverse.getNumIdenticals();
				
				/*System.out.println(leftSideForward);
				System.out.println("totalFwdFlankLength : "+totalForwardFlanksLength);
				System.out.println("totalFwdFlankIdenticals : "+totalForwardFlanksIdenticals);
				
				System.out.println(rightSideReverse);
				System.out.println("totalRevFlankLength : "+totalReverseFlanksLength);
				System.out.println("totalRevFlankIdenticals : "+totalReverseFlanksIdenticals);*/
			}
			else{
				if (rLeftFlankLen<flankSizeLimit && rRightFlankLen>=flankSizeLimit){
					//System.out.println("left is small, right is ok");
					
					totalForwardFlanksLength = rightSideForward.getLength();
					totalForwardFlanksIdenticals = rightSideForward.getNumIdenticals();
					
					totalReverseFlanksLength = leftSideReverse.getLength();
					totalReverseFlanksIdenticals = leftSideReverse.getNumIdenticals();
					
					/*System.out.println(rightSideForward);
					System.out.println("totalFwdFlankLength : "+totalForwardFlanksLength);
					System.out.println("totalFwdFlankIdenticals : "+totalForwardFlanksIdenticals);
					
					System.out.println(leftSideReverse);
					System.out.println("totalRevFlankLength : "+totalReverseFlanksLength);
					System.out.println("totalRevFlankIdenticals : "+totalReverseFlanksIdenticals);*/
					
				}
				else{
					if(rLeftFlankLen<flankSizeLimit && rRightFlankLen<flankSizeLimit){
						//System.out.println("no flanks");
						noFlanks = true;
					}
				}
			}
		}
	
		int totalPairedReadLength = pairedEnd.length();
		int totalPairedEndLeftIdenticals = pairedLeft.getNumIdenticals();
		int totalPairedEndRightIdenticals = pairedRight.getNumIdenticals();
		
		//System.out.println("forward paired end:\n"+pairedRight+(double)totalPairedEndRightIdenticals/totalPairedReadLength);
		
		//System.out.println("reverse paired end:\n"+pairedLeft+(double)totalPairedEndLeftIdenticals/totalPairedReadLength);

		double forwardScoreMin;
		double reverseScoreMin;
		if (noFlanks==false){
			forwardScoreMin = Math.min((double)totalForwardFlanksIdenticals/totalForwardFlanksLength,(double)totalPairedEndRightIdenticals/totalPairedReadLength);
			reverseScoreMin = Math.min((double)totalReverseFlanksIdenticals/totalReverseFlanksLength ,(double)totalPairedEndLeftIdenticals/totalPairedReadLength);
		}else{
			//System.out.println("aligning only paired end sionce the flanks are too short");
			forwardScoreMin=(double)totalPairedEndRightIdenticals/totalPairedReadLength;
			reverseScoreMin=(double)totalPairedEndLeftIdenticals/totalPairedReadLength;
		}
		
		//System.out.println();
		//System.out.println("FINISHED TO ALIGN:");

		if (forwardScoreMin>reverseScoreMin){
			//System.out.println("Forward wins, the minimal score was "+forwardScoreMin);
			return new AlignmentsData(forwardScoreMin, id, readLeftFlank, readRightFlank, readRepetitiveRegion, leftSideForward, rightSideForward, pairedRight);
		}
		else{
			//System.out.println("Reverse wins, the minimal score was "+reverseScoreMin);
			return new AlignmentsData(reverseScoreMin, id, reverseReadRightSeq, reverseReadLeftSeq, reverseReadRepetitiveSeq, leftSideReverse, rightSideReverse, pairedLeft);
		}
	}	
}
