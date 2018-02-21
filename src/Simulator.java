/*
 * This class is our simulator for reads.
 * It takes in a fasta dna file and writes a txt file containing the simulated reads.
 * Author: Yanhan Lyu, Weijia Ma */

import java.io.*;
import java.util.*;

public class Simulator{
    private String dnaFile = "";
    private int coverage = 0;
    private int readLength = 0;
    private double errorRate = 0;

    public Simulator(String dnaFile, int coverage, int readLength, double errorRate){
        this.dnaFile = dnaFile;
        this.coverage = coverage;
        this.readLength = readLength;
        this.errorRate = errorRate;
    }

    /*
     Parse the dna fasta file into a string.
     */
    public String dnaParser(){
        String fileName = this.dnaFile;
        String returnString = "";
        try {
            File file = new File(fileName);
            FileReader fileReader = new FileReader(file);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            StringBuffer stringBuffer = new StringBuffer();
            String line;
            while ((line = bufferedReader.readLine()) != null) {
                if (line.contains(">")){
                    continue;
                }
                stringBuffer.append(line.toUpperCase());

            }
            bufferedReader.close();
            fileReader.close();
            returnString = stringBuffer.toString();
            return stringBuffer.toString();
        } catch (IOException error) {
            error.printStackTrace();
        }
        return returnString;
    }

    /*
     Introduce error to each position of a read
     according to the error rate from input parameter.
     Return a string representing the read with error.
     */
    private String errorHandler(String read){
        String newRead = "";
        for (int i = 0; i < read.length(); i++){
            double probThreshold = Math.random();
            if (probThreshold <= this.errorRate){
                Random random = new Random();
                int randomIndex = random.nextInt(3);
                //System.out.println(randomIndex);
                switch (read.charAt(i))
                {
                    case 'A':
                        char[] errorSetA = {'T','C','G'};
                        newRead = newRead+errorSetA[randomIndex];
                        break;
                    case 'C':
                        char[] errorSetC = {'T','A','G'};
                        newRead = newRead+ errorSetC[randomIndex];
                        break;
                    case 'T':
                        char[] errorSetT = {'C','A','G'};
                        newRead = newRead+ errorSetT[randomIndex];
                        break;
                    case 'G':
                        char[] errorSetG = {'C','A','T'};
                        newRead = newRead+ errorSetG[randomIndex];
                        break;
                    default:
                        break;
                }

            } else {
                newRead = newRead + read.charAt(i);
            }
        }
        return newRead;
    }

    /*
     get the length of dna
     */
    private int getDnaLength(String dna){
        return dna.length();
    }

    /*
     Helper method to calculate the N
     */
    private double getN(int dnaLength, int coverage, int readLength){
        return coverage*dnaLength/readLength;
    }

    /*
     Write the simulated strings of read into a txt file
     */
    private void writeToFile(String reads){
        try {
            PrintWriter writer = new PrintWriter("reads.txt", "UTF-8");
            writer.print(reads);
            writer.close();
        }catch (IOException error) {
            error.printStackTrace();
        }
    }

    /*
     Simulate the reads and output the txt file.
     */
    public void simulation(){
        String dna = dnaParser();
        int dnaLength = getDnaLength(dna);
        double numReads = getN(dnaLength, this.coverage, this.readLength);
        int count = 0;
        String reads = "";
        while (count < numReads){
            Random random = new Random();
            int randomIndex = random.nextInt(dnaLength);
            String subdna = "";
            if (randomIndex + this.readLength < dnaLength) {
                subdna = dna.substring(randomIndex, randomIndex + this.readLength);
            } else {
                subdna = dna.substring(randomIndex);
            }
            String finalRead = errorHandler(subdna);
            reads += finalRead + "\n";
            count ++;
        }
        writeToFile(reads);
    }

    public static void main(String[] args) {
        if (args.length == 4){
            String dnaFile = args[0];
            int coverage = Integer.parseInt(args[1]);
            int readLength = Integer.parseInt(args[2]);
            double errorRate = Double.parseDouble(args[3]);
            Simulator simulator = new Simulator(dnaFile, coverage, readLength, errorRate);
            simulator.simulation();
        }
    }

}

