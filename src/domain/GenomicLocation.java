package domain;

public class GenomicLocation {
	private int chrIndex;
	private int position;
	public GenomicLocation(int chrIndex, int position){
		this.chrIndex = chrIndex;
		this.position = position;
	}
	
	public String toString(){
		return "Temp.rep. "+"\t"+chrIndex+":"+position;
		
	}
	

}
