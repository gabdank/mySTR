package domain;

public class KmerRepeatUnitPair {
	private final int kmerID;
	private final int repeatUnitID;
	
	public int getKmerID(){
		return kmerID;				
	}
	public int getRepeatUnitID(){
		return repeatUnitID;
	}
	
	public KmerRepeatUnitPair(int kmer,int unit){
		kmerID = kmer;
		repeatUnitID = unit;
	}
	
	public boolean equals(Object o){
		int k = ((KmerRepeatUnitPair)o).getKmerID();
		int r = ((KmerRepeatUnitPair)o).getRepeatUnitID();
		return ((k==kmerID) && (repeatUnitID==r));
	}
	public int hashCode(){
		return  Long.valueOf(kmerID * 31 + repeatUnitID).hashCode();
	}
	
	public String toString(){
		return "("+kmerID+":"+repeatUnitID+")";
	}
	
}
