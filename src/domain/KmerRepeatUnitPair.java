package domain;

public class KmerRepeatUnitPair {
	private int kmerID;
	private int repeatUnitID;
	
	public KmerRepeatUnitPair(int kmer,int unit){
		kmerID = kmer;
		repeatUnitID = unit;
	}
	
	public boolean equals(int kmer, int repUnit){
		return kmer==kmerID && repeatUnitID==repUnit;
	}
	
	public int hashCode(){
		return  Long.valueOf(kmerID * 31 + repeatUnitID).hashCode();
	}
	
}
