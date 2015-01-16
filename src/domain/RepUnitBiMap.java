package domain;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import utils.SequenceUtils;


public class RepUnitBiMap {
	private static BiMap<String, Integer> unitID = HashBiMap.create();
	private static BiMap<Integer, String> idUnit = HashBiMap.create();
	
	private static final String s = "AGCT";
	private static int counter = 0;	
	
	public RepUnitBiMap(int maxLengthOfRepeat){
		initiateTable(maxLengthOfRepeat);
		idUnit = unitID.inverse();
	}
	
	/*public static void main(String[] args){
		RepUnitBiMap x = new RepUnitBiMap(3);
		System.out.println(unitID.keySet());
		System.out.println(idUnit.keySet());
	}*/
	
	private void initiateTable(int maxLen){		
		for (int i=2;i<=maxLen;i++){
			permute(i,"",unitID);
		}
				
	}
	
	private void permute(int level, String prefix, BiMap<String, Integer> unitID ){
		if (level == 0){			
			String representative = SequenceUtils.getRepresentative(prefix);
			if (! unitID.containsKey(representative)){
				unitID.put(representative,new Integer(counter));		
				counter++;
			}
			return;
		}
		for (int i = 0; i < s.length(); i++)
			permute(level - 1, prefix + s.charAt(i),unitID);
	}
	
	public int getIntForString(String key){
		if (!unitID.containsKey(key)){
			throw new IllegalArgumentException("The map doesn't contain the k-mer specified ("+key+")...");
		}
		return unitID.get(key);
	}
	
	public String getStringForInt(int key){
		return idUnit.get(key);
	}
	
}
