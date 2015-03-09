package domain;

public class KeyPair {
	private int chrIndex;
	private int chrPosition;
	
	public KeyPair(int ch,int p){
		chrIndex = ch;
		chrPosition = p;
	}

	public boolean equals(Object o){
		return ((KeyPair)o).chrIndex==chrIndex && ((KeyPair)o).chrPosition==chrPosition; 
	}
	
	public int hashCode(){
		int hash = 1;
        hash = hash * 17 + chrIndex;
        hash = hash * 31 + chrPosition;        
        return hash;
		
	}
	
	public int getChrIndex() {
		return chrIndex;
	}

	public void setChrIndex(int chrIndex) {
		this.chrIndex = chrIndex;
	}

	public int getChrPosition() {
		return chrPosition;
	}

	public void setChrPosition(int chrPosition) {
		this.chrPosition = chrPosition;
	}
	
	public String toString(){
		return chrIndex+":"+chrPosition;
	}
	
}
