package domain;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;

import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

public class kmerInformator {
	private KmerMap kmerMap;
	private RepUnitBiMap repeatMap;
	
	private static final int KMER_NON_UNIQUENESS_THRESHOLD = 2;

	private Chromosome2IndexBiMap chromoMap;	
	private int numberOfThreads=-1;

	private TIntObjectMap<ArrayList<GenomicLocation>> repUnitIndex2ListOfLocations;

	private HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>> immidiateFlanksKmerRepUnit2SetOfLocations;
	private HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>> longFlanksKmerRepUnit2SetOfLocations;
	
	private TIntObjectMap<TIntIntMap> repeatUnit2ValidShortKmers;
	private TIntObjectMap<TIntIntMap> repeatUnit2NotValidShortKmers;
	
	private TIntObjectMap<TIntIntMap> repeatUnit2ValidLongKmers;
	private TIntObjectMap<TIntIntMap> repeatUnit2NotValidLongKmers;
	
	

	public kmerInformator(int kLength, String chromosomeListFileName, String trfBedFileName, int numOfThreads) throws IOException{
		numberOfThreads = numOfThreads;
		System.out.println("Beginning generation of the k-mer map (String->int)");
		kmerMap = new KmerMap(kLength);
		System.out.println("Finished generation of the k-mer map (String->int)");

		System.out.println("Beginning generation of the repeat representative map (String->int)");
		repeatMap = new RepUnitBiMap(8);
		System.out.println("Finished generation of the repeat representative map (String->int)");

		System.out.println("Beginning generation of the Chromosome representative map (String->int)");
		chromoMap = new Chromosome2IndexBiMap(chromosomeListFileName);
		System.out.println("Finished generation of the Chromosome representative map (String->int)");

		repUnitIndex2ListOfLocations = new TIntObjectHashMap<ArrayList<GenomicLocation>>(10000);

		immidiateFlanksKmerRepUnit2SetOfLocations = new HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>>(17000000);
		longFlanksKmerRepUnit2SetOfLocations = new HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>>(117000000);
		// helpful maps in order to keep track of un-uniquely mapping k-mers
		
		repeatUnit2ValidLongKmers = new TIntObjectHashMap<TIntIntMap>();
		repeatUnit2NotValidLongKmers = new TIntObjectHashMap<TIntIntMap>();
		
		repeatUnit2ValidShortKmers = new TIntObjectHashMap<TIntIntMap>();
		repeatUnit2NotValidShortKmers = new TIntObjectHashMap<TIntIntMap>();

		System.out.println("STARTING TO PARSE THE BED FILE");
		processTRFBedFile(trfBedFileName);	
		System.out.println("FINISHED TO PARSE THE BED FILE");

	}

	/**
	 * Read trf.bed file that was crerated using different Python scripts
	 * The file contains information about VALID short tandem repeat sites
	 * Each site was validated earlier!!
	 * Anyway - each line in the trf.bed file is a short tandem repeat site.
	 * Each line contains the following information:
	 * 1. representetive repeat unit STRING
	 * 2. chromosome index
	 * 3. start (position in the chromosome)
	 * 4. repeatLength
	 * 5. sequence
	 * 6. startInside
	 * 7. List of all the k-mers coming from the imiidiate flanks (30)
	 * 8. $ - separator
	 * 9. List of all the k-mers coming from the long flanks (600?500 bp long)
	 */
	private void processTRFBedFile(String trfBedFileName) throws IOException{
		BufferedReader  br = new BufferedReader(new FileReader(trfBedFileName));
		String line;
		int mone =0;

		long startTime = System.currentTimeMillis();
		long endTime = System.currentTimeMillis();

		while ((line=br.readLine())!=null){
			mone++;
			if (mone%5==0) {
				endTime = System.currentTimeMillis();
				System.out.println("PROCESSED "+mone +" lines in BED, it took "+((endTime-startTime)/1000) +" seconds" );
				startTime = endTime;
				break;
			}
			


			StringTokenizer t = new StringTokenizer(line);
			if (t.countTokens()>4){
				String representativeUnit = t.nextToken();
				int chromosomeIndex = Integer.parseInt(t.nextToken());
				int position = Integer.parseInt(t.nextToken());
				int repeatLength = Integer.parseInt(t.nextToken());

				String fullSeq = t.nextToken();				
				int repeatStart = Integer.parseInt(t.nextToken());

				GenomicLocation gl = new GenomicLocation(chromosomeIndex,position);
				//chromosomeIndex,position,fullSeq, representativeUnit.length(), repeatStart, repeatLength, numberOfThreads);				
				
				//System.out.println("Generated genomic location:");
				//System.out.println(gl);

				addRepeat2LocationMapping(representativeUnit,gl);
				System.out.println(representativeUnit);
				System.out.println(repUnitIndex2ListOfLocations.get(repeatMap.getIntForString("AAGCCT")));
				
				//System.out.println(repeatMap.);
				int repeatUnitIndex = repeatMap.getIntForString(representativeUnit);

				// In order to not take into consederation k-mers that are repetitive and map all over the place - we are requiring some level of 
				// uniqueness. It means that currently any pair of ([k-mer],[representative repeat unit index]) that maps to over 5 
				// genomic locations will be discarded - as not helpful for the mapping purposes.
				// in order to be able to count the references of each pair we have to create a data structure that will keep them
				
				
				TIntIntMap validShortKmersMap = repeatUnit2ValidShortKmers.get(repeatUnitIndex);
				TIntIntMap notValidShortKmersMap = repeatUnit2NotValidShortKmers.get(repeatUnitIndex);

				
				if (validShortKmersMap == null){
					validShortKmersMap = new TIntIntHashMap();
					repeatUnit2ValidShortKmers.put(repeatUnitIndex,validShortKmersMap);
				}
				if (notValidShortKmersMap == null){
					notValidShortKmersMap = new TIntIntHashMap();
					repeatUnit2NotValidShortKmers.put(repeatUnitIndex,notValidShortKmersMap);
				}
				
				System.out.println("validShortKmers:"+validShortKmersMap);
				System.out.println("notValidShortKmers:"+notValidShortKmersMap);
				
				String nextToken = t.nextToken();
				while (!nextToken.equals("$")){
					int kmerIndexValue = Integer.parseInt(nextToken);
					
					if (!notValidShortKmersMap.containsKey(kmerIndexValue)){ // only in case the kmer is potentialy valid
						KmerRepeatUnitPair mappingPair = new KmerRepeatUnitPair(kmerIndexValue,repeatUnitIndex);
						if (validShortKmersMap.containsKey(kmerIndexValue)){ // not the first time encountering kmer
							
							int numOfOccurrences = validShortKmersMap.get(kmerIndexValue);
							if (numOfOccurrences>KMER_NON_UNIQUENESS_THRESHOLD){
								notValidShortKmersMap.put(kmerIndexValue, 1);
								validShortKmersMap.remove(kmerIndexValue);
								// remove references to kmer from the locations?
								System.out.println("REMOVING FROM SHORT!!!");
								removeKmerUnitPairFromMapping(mappingPair,immidiateFlanksKmerRepUnit2SetOfLocations );
								//removeKmerUnit2ImmidiateLocationMapping(mappingPair);


							}
							else{
								// add reference to the kmer to the locations...
								validShortKmersMap.increment(kmerIndexValue);
								//addKmer2ShortLocationMapping(kmerIndexValue, gl);
								addKmerUnitPairFromMapping(mappingPair,gl,immidiateFlanksKmerRepUnit2SetOfLocations );

							}
						}else{ // first time encountering the kmer
							// add reference to the kmer to the locations...
							validShortKmersMap.put(kmerIndexValue, 1);
							addKmerUnitPairFromMapping(mappingPair,gl, immidiateFlanksKmerRepUnit2SetOfLocations );						}
					}
					nextToken = t.nextToken();
				}

				System.out.println("Finished to do the short flanks kmers");
				System.out.println("Key SET: "+immidiateFlanksKmerRepUnit2SetOfLocations.keySet());
				for (KmerRepeatUnitPair p : immidiateFlanksKmerRepUnit2SetOfLocations.keySet()){
					System.out.println(p+" : "+immidiateFlanksKmerRepUnit2SetOfLocations.get(p));
					
				}
				//System.out.println();
				//System.out.println("CHECKING:"+immidiateFlanksKmerRepUnit2SetOfLocations.get(new KmerRepeatUnitPair(23, 242)));
				//////////////////////////////////////////
				TIntIntMap validLongKmersMap = repeatUnit2ValidLongKmers.get(repeatUnitIndex);
				TIntIntMap notValidLongKmersMap = repeatUnit2NotValidLongKmers.get(repeatUnitIndex);
				
				
				if (validLongKmersMap == null){
					validLongKmersMap = new TIntIntHashMap();
					repeatUnit2ValidLongKmers.put(repeatUnitIndex,validLongKmersMap);
				}
				if (notValidLongKmersMap == null){
					notValidLongKmersMap = new TIntIntHashMap();
					repeatUnit2NotValidLongKmers.put(repeatUnitIndex,notValidLongKmersMap);
				}

				while (t.hasMoreTokens()){
					int kmerIndexValue = Integer.parseInt(t.nextToken());
					
					if (!notValidLongKmersMap.containsKey(kmerIndexValue)){ // only in case the kmer is potentialy valid
						KmerRepeatUnitPair mappingPair = new KmerRepeatUnitPair(kmerIndexValue,repeatUnitIndex);

						if (validLongKmersMap.containsKey(kmerIndexValue)){ // not the first time encountering kmer
							int numOfOccurrences = validLongKmersMap.get(kmerIndexValue);
							if (numOfOccurrences>KMER_NON_UNIQUENESS_THRESHOLD){
								System.out.println("REMOVING FROM LONG!!!");
								notValidLongKmersMap.put(kmerIndexValue, 1);
								validLongKmersMap.remove(kmerIndexValue);
								
								removeKmerUnitPairFromMapping(mappingPair,longFlanksKmerRepUnit2SetOfLocations );

							}
							else{
								// add reference to the kmer to the locations...
								validLongKmersMap.increment(kmerIndexValue);
								addKmerUnitPairFromMapping(mappingPair,gl, longFlanksKmerRepUnit2SetOfLocations );	

							}
						}else{ // first time encountering the kmer
							// add reference to the kmer to the locations...
							validLongKmersMap.put(kmerIndexValue, 1);
							addKmerUnitPairFromMapping(mappingPair,gl, longFlanksKmerRepUnit2SetOfLocations );	
						}
					}
				}
				System.out.println("Finished to do the long flanks kmers");
				System.out.println("Size of the long list:"+longFlanksKmerRepUnit2SetOfLocations.keySet().size());
				System.out.println(longFlanksKmerRepUnit2SetOfLocations.keySet());
				System.out.println();
			}
		}
	}

	private void addRepeat2LocationMapping(String repeatUnitRepresentative,GenomicLocation loc){
		addIndexLocationToMap(repeatMap.getIntForString(repeatUnitRepresentative), loc, repUnitIndex2ListOfLocations);
	}

	private void addIndexLocationToMap(int index, GenomicLocation l, TIntObjectMap<ArrayList<GenomicLocation>> m){
		if (m.containsKey(index)){
			ArrayList<GenomicLocation> l1 = (ArrayList<GenomicLocation>) m.get(index);
			l1.add(l);

		}else{
			ArrayList<GenomicLocation> l1 = new ArrayList<GenomicLocation>();
			l1.add(l);
			m.put(index, l1);
		}
	}
	
	private void removeKmerUnitPairFromMapping(KmerRepeatUnitPair mappingPair, HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>> mapping){
		mapping.remove(mappingPair);
	}
	
	private void addKmerUnitPairFromMapping(KmerRepeatUnitPair mappingPair,GenomicLocation loc, HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>>  mapping){
		if (mapping.containsKey(mappingPair)){
			HashSet<GenomicLocation> l1 = mapping.get(mappingPair);
			if (!l1.contains(loc)){
				l1.add(loc);
			}
		}
		else{
			HashSet<GenomicLocation> l1 = new HashSet<GenomicLocation>();
			l1.add(loc);
			mapping.put(mappingPair, l1);
		}
			
	}
}



