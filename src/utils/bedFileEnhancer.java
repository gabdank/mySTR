package utils;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.hash.TIntIntHashMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import java.util.StringTokenizer;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import domain.KmerMap;

public class bedFileEnhancer {
	public static void main(String[] args) throws IOException{
		int kmerSize=12;
		System.out.println("STARTEDED to create 12-mer table");
		KmerMap kmersTable = new KmerMap(kmerSize);
		System.out.println("FINISHED to create 12-mer table");

		/*-javaagent:/home/gabdank/Programs/SizeOf_0_2_1/SizeOf.jar
		 * SizeOf.skipStaticField(true); //java.sizeOf will not compute static fields
		SizeOf.skipFinalField(true); //java.sizeOf will not compute final fields
		SizeOf.skipFlyweightObject(true); //java.sizeOf will not compute well-known flyweight objects
		System.out.println(SizeOf.humanReadable(SizeOf.deepSizeOf((kmersTable)))); //this will print the object size in bytes
		 */

		BufferedReader chrReader = new BufferedReader(new FileReader("/home/gabdank/Documents/HumanSTR/chromosomes.list"));



		BiMap<String, Integer> chromoID = HashBiMap.create();


		String word = "";
		int mone =0;
		while ((word=chrReader.readLine())!=null){
			StringTokenizer t = new StringTokenizer(word);
			String toPut = t.nextToken();

			//chromosomes.put(toPut,new Integer(mone));
			chromoID.put(toPut, new Integer(mone));

			//System.out.println();
			mone++;
		}
		chrReader.close();

		
		

		BufferedReader  br = new BufferedReader(new FileReader("/home/gabdank/Documents/HumanSTR/hg19.two.filtered.bed"));
		//BufferedReader  br = new BufferedReader(new FileReader("trfIndexNEW1.bed"));
		String line;
		mone =0;

		System.out.println("STARTED READING BED FILE");


		FileWriter fw1 = new FileWriter(new File("/home/gabdank/Documents/HumanSTR/hg19.filtered.indexed.bed"));
		BufferedWriter bw1 = new BufferedWriter(fw1);


		while ((line=br.readLine())!=null){
			if (mone%100000==0) 
				System.out.println("PROCESSED "+mone +" lines in BED");
			mone++;
			/*if (mone==10){
				break;
			}*/
			StringTokenizer t = new StringTokenizer(line);
			int numberOfTokens = t.countTokens();
			//System.out.println("Number of tokens: "+ numberOfTokens);
			if (numberOfTokens == 7){
				String chromosome = t.nextToken();
				int start = Integer.parseInt(t.nextToken());
				int end = Integer.parseInt(t.nextToken());

				String unit=t.nextToken();
				String sequence=t.nextToken();
				int repLength = Integer.parseInt(t.nextToken());
				int startInside = Integer.parseInt(t.nextToken());
				
				String representativeUnit = SequenceUtils.getRepresentative(unit);

				int chromosomeIndex = chromoID.get(chromosome).intValue();
				int repeatLength = end-start+1;

				/**
				 * We have to create kmers for the left and right short flanks (about 20-30 nucleotides) to 
				 * allow relatively fast way of filtering. The problem is obviously with reads that are
				 * almost fully repetitive. In these cases the paired end is going to be useful. 
				 * So I am creating two groups of kmer - close to the repeat and paired end k-mers (that are 100-200bp away from the repeat.
				 * 
				 * sequence contains left flank (that could be relatively short, and same true to the right flank.
				 * 
				 **/
				String leftFlank = sequence.substring(0,startInside-1);
				String rightFlank = sequence.substring(startInside-1+repLength);
				
				String shortLeft = "";
				String shortRight = "";
						
				if (leftFlank.length()<30){
					shortLeft = leftFlank;
				}else{
					shortLeft = leftFlank.substring(leftFlank.length()-30);
				}
				
				if (rightFlank.length()<30){
					shortRight = rightFlank;
				}else{
					shortRight = rightFlank.substring(0,30);
				}		
				
				TIntIntHashMap shortFlanksMap = new TIntIntHashMap(1000);
				
				for (int i=0;i<shortLeft.length()-kmerSize+1;i++){
					String twelveMer = shortLeft.substring(i,i+kmerSize);				
					if (!twelveMer.contains("N")){			
						shortFlanksMap.putIfAbsent(kmersTable.getIntForString(twelveMer),1);					
						String revComp = SequenceUtils.getReverseComplementary(twelveMer);
						shortFlanksMap.putIfAbsent(kmersTable.getIntForString(revComp),1);
					}
				}
				for (int i=0;i<shortRight.length()-kmerSize+1;i++){
					String twelveMer = shortRight.substring(i,i+kmerSize);				
					if (!twelveMer.contains("N")){			
						shortFlanksMap.putIfAbsent(kmersTable.getIntForString(twelveMer),1);					
						String revComp = SequenceUtils.getReverseComplementary(twelveMer);						
						shortFlanksMap.putIfAbsent(kmersTable.getIntForString(revComp),1);
					}
				}

				TIntIntHashMap longFlanksMap = new TIntIntHashMap(1000);
				for (int i=0;i<leftFlank.length()-kmerSize+1;i++){
					String twelveMer = leftFlank.substring(i,i+kmerSize);				
					if (!twelveMer.contains("N")){			
						longFlanksMap.putIfAbsent(kmersTable.getIntForString(twelveMer),1);					
						String revComp = SequenceUtils.getReverseComplementary(twelveMer);
						longFlanksMap.putIfAbsent(kmersTable.getIntForString(revComp),1);
					}
				}
				for (int i=0;i<rightFlank.length()-kmerSize+1;i++){
					String twelveMer = rightFlank.substring(i,i+kmerSize);				
					if (!twelveMer.contains("N")){			
						longFlanksMap.putIfAbsent(kmersTable.getIntForString(twelveMer),1);					
						String revComp = SequenceUtils.getReverseComplementary(twelveMer);
						longFlanksMap.putIfAbsent(kmersTable.getIntForString(revComp),1);
					}
				}
				
				TIntIterator iter =  shortFlanksMap.keySet().iterator();

				bw1.write(representativeUnit+"\t"+chromosomeIndex+"\t"+start+"\t"+repeatLength);
				bw1.write("\t"+sequence+"\t"+startInside);
				while (iter.hasNext()){
					bw1.write("\t"+iter.next());
				}
				bw1.write("\t$\t");
				
				iter =  longFlanksMap.keySet().iterator();
				
				while (iter.hasNext()){
					bw1.write("\t"+iter.next());
				}
				bw1.write("\n");	
				
			}
			else{
				System.out.println(line);
			}
		}	
		System.out.println("Mone " + mone);
		System.out.println("FINISHED READING BED FILE");
		bw1.close();
		fw1.close();
	}
}
