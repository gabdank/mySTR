package test;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

import domain.kmerInformator;

public class KmerRepeatUnitPairTest {

	@Test
	public void test() throws IOException {
		kmerInformator testInformator = new kmerInformator(12, "/home/gabdank/Documents/STR_Attempt/chromosomes.list", "/home/gabdank/Documents/STR_Attempt/temp.indexed.bed", 5);
	}

}
