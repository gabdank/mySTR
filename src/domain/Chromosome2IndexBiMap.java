package domain;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

public class Chromosome2IndexBiMap {
	private BiMap<String, Integer> chromoID = HashBiMap.create();
	private BiMap<Integer, String> idChromo = HashBiMap.create();
	
	public Chromosome2IndexBiMap(String chromosomeListFileName){
		try {
			BufferedReader chrReader = new BufferedReader(new FileReader(chromosomeListFileName));
			String word = "";
			int mone =0;
			while ((word=chrReader.readLine())!=null){
				StringTokenizer t = new StringTokenizer(word);
				String toPut = t.nextToken();
				chromoID.put(toPut, new Integer(mone));
				mone++;
			}
			chrReader.close();
			idChromo = chromoID.inverse();
			
			
		} catch (FileNotFoundException e) {
			System.out.println("No File containing chromosomes list was found...");
			e.printStackTrace();
		} catch (IOException e) {
			System.out.println("Problem reading the file containing chromosomes list...");
			e.printStackTrace();
		}

	}
	
	public int getChromosomeIndex(String c){
		return chromoID.get(c);
	}
	
	public String getChroosomeName(int index){
		return idChromo.get(new Integer(index));
	}
}
