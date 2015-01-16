package domain;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import gnu.trove.map.TIntObjectMap;

public class kmerInformator {
	private KmerMap kmerMap;
	private RepUnitBiMap repeatMap;

	private Chromosome2IndexBiMap chromoMap;	
	private final int numberOfThreads=-1;
	
	private TIntObjectMap<ArrayList<GenomicLocation>> repUnitIndex2ListOfLocations;
	private HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>> immidiateFlanksKmerRepUnit2SetOfLocations;
	private HashMap<KmerRepeatUnitPair,HashSet<GenomicLocation>> longFlanksKmerRepUnit2SetOfLocations;
	
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

		shortKmerIndex2ListOfLocations = new TIntObjectHashMap<HashSet<GenomicLocation>>(17000000);
		
		longKmerIndex2ListOfLocations = new TIntObjectHashMap<HashSet<GenomicLocation>>(17000000);
		
		repUnitIndex2ListOfLocations = new TIntObjectHashMap<ArrayList<GenomicLocation>>(10000);

		//uniqueKmers = new TIntIntHashMap();
		//notUniqueKmers = new TIntIntHashMap();

		repeatUnit2ValidLongKmers = new TIntObjectHashMap<TIntIntMap>();
		repeatUnit2NotValidLongKmers = new TIntObjectHashMap<TIntIntMap>();
		
		repeatUnit2ValidShortKmers = new TIntObjectHashMap<TIntIntMap>();
		repeatUnit2NotValidShortKmers = new TIntObjectHashMap<TIntIntMap>();
		
		System.out.println("STARTING TO PARSE THE BED FILE");
		processTRFBedFile(trfBedFileName);	
		System.out.println("FINISHED TO PARSE THE BED FILE");

		int[] kmers = longKmerIndex2ListOfLocations.keys();
		int maxValue = 0;
		for (int i=0;i<kmers.length;i++){
			if (longKmerIndex2ListOfLocations.get(kmers[i]).size()>maxValue){
				maxValue = longKmerIndex2ListOfLocations.get(kmers[i]).size();
			}			
		}	
		
		int[] occurrences = new int[maxValue+10];
		for (int i=0;i<kmers.length;i++){
			if (longKmerIndex2ListOfLocations.get(kmers[i]).size()>0){
				occurrences[longKmerIndex2ListOfLocations.get(kmers[i]).size()]++;
			}					
		}
		
	}
	
	
}



