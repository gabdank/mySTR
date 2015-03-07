package test;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import org.junit.Test;

import domain.GenomicLocation;
import domain.KmerRepeatUnitPair;
import domain.kmerInformator;
import domain.multiThreading.processingThread;

public class ThreadTest {

	@Test
	public void test() throws IOException {
		kmerInformator testInformator = new kmerInformator(12, "/home/gabdank/Documents/STR_Attempt/Simulation/chromosomes.list", "/home/gabdank/Documents/STR_Attempt/Simulation/one.indexed.bed", 13);
		HashMap<KmerRepeatUnitPair, HashSet<GenomicLocation>> z = testInformator.getImmidiateFlanksKmerRepUnit2SetOfLocations();
		/*for (KmerRepeatUnitPair p : z.keySet()){
			System.out.println(p+"\t"+z.get(p));
		}*/
		
		processingThread t = new processingThread(5,"/home/gabdank/Documents/STR_Attempt/Simulation/case7_1.fastq","/home/gabdank/Documents/STR_Attempt/Simulation/case7_2.fastq", testInformator,"/home/gabdank/Documents/STR_Attempt/Simulation/out_") ;
		t.run();
	}
}
