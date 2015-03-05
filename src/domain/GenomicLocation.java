package domain;

public class GenomicLocation {
	private int chrIndex;
	private int position;
	
	private final String fullSequence;
	private final int repeatUnitLength;
	private final int repetitiveSequenceStartPosition;	
	private final int repetitiveSequenceLength;
	
	private int[] multyThreadCounters;

	public GenomicLocation(int chrIndex, int position, String fullSeq, int repUnitLength, int repeatStartPosition, int repeatLength, int numberOfThreads){
		this.chrIndex = chrIndex;
		this.position = position;
		fullSequence=fullSeq;
		repeatUnitLength = repUnitLength;
		repetitiveSequenceStartPosition = repeatStartPosition;
		repetitiveSequenceLength = repeatLength;
		multyThreadCounters = new int[numberOfThreads];

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

	public int[] getMultyThreadCounters() {
		return multyThreadCounters;
	}

	public void setMultyThreadCounters(int[] multyThreadCounters) {
		this.multyThreadCounters = multyThreadCounters;
	}

	public int getChromosome() {
		return chrIndex;
	}

	public int getStartPosition() {
		return position;
	}

	public String getFullSequence() {
		return fullSequence;
	}

	public int getRepeatUnitLength() {
		return repeatUnitLength;
	}

	public int getRepetitiveSequenceStartPosition() {
		return repetitiveSequenceStartPosition;
	}

	public int getRepetitiveSequenceLength() {
		return repetitiveSequenceLength;
	}
	public String toString(){
		return "GenomicLocation:chrIndex_"+chrIndex+":position_"+position+":"+
	fullSequence.substring(repetitiveSequenceStartPosition,repetitiveSequenceStartPosition+repeatUnitLength);
	}
	
	 
	
	public String reportTheGenomicLocation(){
		String toReturn = "";
		toReturn = "repetitiveUnit:"+getRepetitiveSection().substring(0, repeatUnitLength)+"$"+"fullSequence:"+fullSequence+"$leftFlank:"+getLeftFlank()+
				"$rightFlank:"+getRightFlank()+"$repetitiveSequence:"+getRepetitiveSection()+
				"$chromosome:"+chrIndex+"$position:"+position+"$repetitiveLength:"+repetitiveSequenceLength+"$";
		
		return toReturn;
	}
	
	
	public String getLeftFlank(){
		return fullSequence.substring(0, repetitiveSequenceStartPosition-1);
	}
	
	public String getRightFlank(){
		return fullSequence.substring(repetitiveSequenceStartPosition-1+repetitiveSequenceLength);
	}
	
	public String getRepetitiveSection(){
		return fullSequence.substring(repetitiveSequenceStartPosition-1,repetitiveSequenceStartPosition-1+repetitiveSequenceLength);
	}

}
