package domain;

import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class KmerMap {
	private static TObjectIntMap<String> kmersMap;
	private static String s = "AGCT";
	private final int k;
	private static int counter = 0;	
	public KmerMap(int kMerLength){
		k=kMerLength;
		initiateTable();
	}	
	private void initiateTable(){
		int quantity = (int) Math.pow(4,k);
		kmersMap = new TObjectIntHashMap<String>(quantity);
		permute(k,"",kmersMap);		
	}
	private void permute(int level, String prefix,  TObjectIntMap<String>  t){
		if (level == 0){			
			t.put(prefix,counter);			
			counter++;
			return;
		}
		for (int i = 0; i < s.length(); i++)
			permute(level - 1, prefix + s.charAt(i),t);
	}
	public int getIntForString(String key){
		if (!kmersMap.containsKey(key)){
			throw new IllegalArgumentException("The map doesn't contain the k-mer specified ("+key+")...");
		}
		return kmersMap.get(key);
	}
}
