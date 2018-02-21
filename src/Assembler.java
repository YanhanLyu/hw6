/*
 * This class is our assembler.
 * It takes in a txt file of reads and writes out a txt file of contigs.
 * Currently, it prints out the de Bruijn edges.
 * Author: Yanhan Lyu, Weijia Ma */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.*;

public class Assembler {
    private String readsFile = "";
    private int k = 0;

    public Assembler(String readsFile, int k){
        this.readsFile = readsFile;
        this.k = k;
    }

    /*
     Parse the txt file of reads into an array list of strings of reads
     */
    public ArrayList<String> readsParser(){
        String fileName = this.readsFile;
        ArrayList<String> readsList = new ArrayList<String>();
        try {
            File file = new File(fileName);
            FileReader fileReader = new FileReader(file);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String read;
            while ((read = bufferedReader.readLine()) != null) {
                if (read.contains(">")){
                    continue;
                }
                readsList.add(read.toUpperCase());
            }
            bufferedReader.close();
            fileReader.close();
            return readsList;
        } catch (IOException error) {
            error.printStackTrace();
        }
        return readsList;
    }

    /*
     Generate an array list of k-mers from the reads list
     */
    private ArrayList<String> formKMers(ArrayList<String> readsList) {
        ArrayList<String> kMers = new ArrayList<String>();
        for (String read : readsList){
            int n = read.length();
            int i = 0;
            // what about reads < k? right now we are ignoring them
            while (i+k-1<n){
                kMers.add(read.substring(i,i+k));
                i++;
            }
        }
        return kMers;
    }


    /*
     Generate an array list of k-mers from the reads list
     */
     private ArrayList<String> formKMinus1Mers(ArrayList<String> kMers) {
        ArrayList<String> KMinus1Mers = new ArrayList<String>();
        for (String kMer : kMers){
            int n = kMer.length();
            int i = 0;
            // what about reads < k? right now we are ignoring them
            KMinus1Mers.add(kMer.substring(0,n-1));
            KMinus1Mers.add(kMer.substring(1,n));
        }
        return KMinus1Mers;
    }

    /*
     Create the de Bruijn edges
     */
    private void assemble() {
        ArrayList<String> readsList = readsParser();
        ArrayList<String> kMers = formKMers(readsList);
        ArrayList<String> kMinus1Mers = formKMinus1Mers(kMers);
        // the de Bruijn edges
        HashMap<Node, ArrayList<Node>> edges = new HashMap<Node, ArrayList<Node>>();
        // dictionary matching a (k-1)mer
        // to its respective node (vertex in the de Bruijn edges)
        HashMap<String, Node> vertices = new HashMap<String, Node>();

        for (int i = 0; i < kMinus1Mers.size()-1; i++){
            String kmerLeft = kMinus1Mers.get(i);
            String kmerRight = kMinus1Mers.get(i+1);
            Node nodeL;
            Node nodeR;
            if (vertices.containsKey(kmerLeft)){
                nodeL = vertices.get(kmerLeft);
            } else {
                nodeL = new Node(kmerLeft);
                vertices.put(kmerLeft,nodeL);
            }

            if (vertices.containsKey(kmerRight)){
                nodeR = vertices.get(kmerRight);
            } else {
                nodeR = new Node(kmerRight);
                vertices.put(kmerRight,nodeR);
            }
            if (edges.containsKey(nodeL)){
                edges.get(nodeL).add(nodeR);
            } else {
                ArrayList<Node> list = new ArrayList<Node>();
                list.add(nodeR);
                edges.put(nodeL, list);
            }
        }
        printEdges(edges);
    }

    /*
     Help us visualize the edges of the de Bruijn edges
     */
    private void printEdges(HashMap<Node, ArrayList<Node>> edges){
        for(Node nodeL:edges.keySet()) {
            for (int i = 0; i < edges.get(nodeL).size(); i ++){
                System.out.print(nodeL.data+" -> ");
                System.out.print(edges.get(nodeL).get(i).data+"\n");
            }
        }
    }

    public static void main(String[] args) {
        if (args.length == 2){
            String readsFile = args[0];
            int k = Integer.parseInt(args[1]);
            Assembler assembler = new Assembler(readsFile, k);
            assembler.assemble();
        }
    }
}
