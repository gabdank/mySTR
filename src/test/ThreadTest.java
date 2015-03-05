package test;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

import domain.kmerInformator;
import domain.multiThreading.processingThread;

public class ThreadTest {

	@Test
	public void test() throws IOException {
		kmerInformator testInformator = new kmerInformator(12, "/home/gabdank/Documents/STR_Attempt/Simulation/chromosomes.list", "/home/gabdank/Documents/STR_Attempt/Simulation/sample.indexed.bed", 13);
		processingThread t = new processingThread(5,"/home/gabdank/Documents/STR_Attempt/Simulation/paired_dat1.fastq","/home/gabdank/Documents/STR_Attempt/Simulation/paired_dat2.fastq", testInformator,"/home/gabdank/Documents/STR_Attempt/Simulation/out_") ;
		t.run();
	}
}
