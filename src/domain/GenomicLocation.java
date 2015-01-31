package domain;

public class GenomicLocation {
	private int chrIndex;
	private int position;
	
	private int[] multyThreadCounters;

	public GenomicLocation(int chrIndex, int position, int numberOfThreads){
		this.chrIndex = chrIndex;
		this.position = position;
		multyThreadCounters = new int[numberOfThreads];

	}
	
	public String toString(){
		return "Temp.rep. "+"\t"+chrIndex+":"+position;
		
	}
	
	public void addOneToCounter(int index){
		//System.out.println("ADDING TO THE "+index + ", while the length of counters is: "+ countersMT.length );
		multyThreadCounters[index]++;		
	}
	
	public void nullifyCounter(int index){
		multyThreadCounters[index] = 0;		
	}
	public int getCounter(int index){
		if (index>=multyThreadCounters.length || index<0){
			System.out.println("Access to negative counter of some thread in GenomicLocation");
			return -1;
		}
		return multyThreadCounters[index];
	
	}
}
