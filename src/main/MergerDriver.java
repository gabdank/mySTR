package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import domain.Filter;
import domain.KeyPair;

public class MergerDriver {


	public static void main(String[] args) throws IOException {
		
		HashMap<KeyPair,ArrayList<String>> alignments = new HashMap<KeyPair,ArrayList<String>>();
		HashMap<KeyPair,String> genomicLocations = new HashMap<KeyPair,String>();
		String direct = "/media/gabdank/Disk3/mySTR/MY14/"; 
		File directory = new File(direct);
		File[] files = directory.listFiles(new Filter("ou"));

		BufferedReader[]  readers = new BufferedReader[files.length];
		for (int i = 0; i<readers.length;i++){
			//String p = .getPath();
			readers[i] = new BufferedReader(new FileReader(files[i]));
		}

		
		String l;
		

		for (int readerId = 0;readerId<readers.length;readerId++){
			int mone = 0;
			int chrInd = 0;
			int posInd = 0;
			while ((l = readers[readerId].readLine())!=null){
				mone++;
				if (mone == 1){
					StringTokenizer tokeniz= new StringTokenizer(l,"$");
					boolean readMore = true;
					while (tokeniz.hasMoreTokens() && readMore){
						// first there is data of the genomic location, and only after that the alignment. 
						// here we will initiate a HashMap of all locations and will "attach" the alignments data to them
						String token = tokeniz.nextToken();
						//System.out.println("TOKEN = "+token);
						StringTokenizer internal = new StringTokenizer(token,":");
						String name = internal.nextToken();
						String val = internal.nextToken();
						if (name.equals("chromosome")){
							chrInd = Integer.parseInt(val);
						}
						if (name.equals("position")){
							posInd = Integer.parseInt(val);
							KeyPair k = new KeyPair(chrInd,posInd);
							
							StringTokenizer important = new StringTokenizer(l,"{");
							String glDATA = important.nextToken();
							
							String DATA = important.nextToken();
							String alignmentData = DATA.substring(0,DATA.length())+"]";

							if (!genomicLocations.containsKey(k)){
								genomicLocations.put(k,glDATA);
							}
							
							if (alignments.containsKey(k)){
								alignments.get(k).add(alignmentData);
							}else{
								ArrayList<String> alis = new ArrayList<String>();
								alis.add(alignmentData);
								alignments.put(k,alis);
							}
							readMore = false;
						}

					}
				}
				if (mone == 2){// alignments data.. may be necessary in the long run
					mone = 0;
				}		

			}
		}



		for (int i = 0; i<readers.length;i++){
			readers[i].close();
		}

		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(direct+"merged.output")));
		for (KeyPair k : alignments.keySet()){
			bw.write(genomicLocations.get(k)+"\n");
			bw.write("{\n");			
			for (int alInd=0;alInd<alignments.get(k).size();alInd++){
				bw.write(alignments.get(k).get(alInd)+"\n");
			}
			
			bw.write("}\n");
		}
		bw.close();

	}

}
