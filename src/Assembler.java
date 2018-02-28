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
    private HashMap<String, ArrayList<String>> formKMers(ArrayList<String> readsList) {
        HashMap<String, ArrayList<String>> readskMers = new HashMap<String, ArrayList<String>>();
        for (String read : readsList){
            int n = read.length();
            int i = 0;
            // what about reads < k? right now we are ignoring them

            while (i+k-1<n){
                if (readskMers.containsKey(read)){
                    readskMers.get(read).add(read.substring(i,i+k));
                } else {
                    ArrayList<String> list = new ArrayList<>();
                    readskMers.put(read, list);
                    list.add(read.substring(i,i+k));
                }
                //System.out.println(read.substring(i,i+k));
                i++;
            }
        }

//        for (Map.Entry<String, ArrayList<String>> entry: readskMers.entrySet()) {
//            String read = entry.getKey();
//            ArrayList<String> kmers = entry.getValue();
//            for (int i = 0; i < kmers.size(); i ++){
//                System.out.println(read+" "+kmers.get(i));
//            }
//        }
        return readskMers;
    }


    /*
     Generate an array list of k-mers from the reads list
     */
    private HashMap<String, ArrayList<String>> formKMinus1Mers(HashMap<String, ArrayList<String>> kMers) {
        HashMap<String, ArrayList<String>>  kMinus1Mers = new HashMap<String, ArrayList<String>>();
        for (Map.Entry<String, ArrayList<String>> entry: kMers.entrySet()) {
            String read = entry.getKey();
            ArrayList<String> kmers = entry.getValue();
            for (String kmer : kmers) {
                int n = kmer.length();
                int i = 0;
                // what about reads < k? right now we are ignoring them
                if (kMinus1Mers.containsKey(read)) {
                    kMinus1Mers.get(read).add(kmer.substring(0, n - 1));
                    kMinus1Mers.get(read).add(kmer.substring(1, n));
                } else {
                    ArrayList<String> list = new ArrayList<>();
                    kMinus1Mers.put(read, list);
                    list.add(kmer.substring(0, n - 1));
                    list.add(kmer.substring(1, n));
                }
            }

        }

//        for (Map.Entry<String, ArrayList<String>> entry: kMinus1Mers.entrySet()) {
//            String read = entry.getKey();
//            ArrayList<String> kmers = entry.getValue();
//            for (int i = 0; i < kmers.size(); i ++){
//                System.out.println(read+" "+kmers.get(i));
//            }
//        }
        return kMinus1Mers;
    }

    /*
     Create the de Bruijn edges
     */
    private void assemble() {
        ArrayList<String> readsList = readsParser();
        HashMap<String, ArrayList<String>> readsKmer = formKMers(readsList);
        HashMap<String, ArrayList<String>> readskMinus1Mers = formKMinus1Mers(readsKmer);
        // the de Bruijn edges
        HashMap<Node, HashSet<Node>> edges = new HashMap<Node, HashSet<Node>>();
        // dictionary matching a (k-1)mer string
        // to its respective node (vertex in the de Bruijn edges)
        HashMap<String, Node> vertices = new HashMap<String, Node>();


        for (Map.Entry<String, ArrayList<String>> entry : readskMinus1Mers.entrySet()) {

            String read = entry.getKey();
            ArrayList<String> kMinus1Mers = entry.getValue();
//            int count = 0;
//            String headkminusmer = "";
//            String endkminusmer = "";

            for (int i = 0; i < kMinus1Mers.size(); i = i + 2) {

//                headkminusmer = kMinus1Mers.get(0);
//                endkminusmer = kMinus1Mers.get(kMinus1Mers.size()-1);
                String kmerLeft = kMinus1Mers.get(i);
                String kmerRight = kMinus1Mers.get(i + 1);
//                System.out.println("??? "+read+"\n");
//                System.out.println("left: "+kmerLeft+"\n");
//                System.out.println("right: "+kmerRight+"\n");
//                System.out.println("begin: "+headkminusmer+"\n");
//                System.out.println("end: "+endkminusmer+"\n");

                //System.out.println("l: "+kmerLeft);
                //System.out.println("r: "+kmerRight);
                Node nodeL;
                Node nodeR;
                if (vertices.containsKey(kmerLeft)) {
                    nodeL = vertices.get(kmerLeft);
                    nodeL.addReads(read);
                } else {
                    nodeL = new Node(kmerLeft,read);
                    vertices.put(kmerLeft, nodeL);
                }

                if (vertices.containsKey(kmerRight)) {
                    nodeR = vertices.get(kmerRight);
                    nodeR.addReads(read);
                } else {
                    nodeR = new Node(kmerRight,read);
                    vertices.put(kmerRight, nodeR);
                }
//                nodeL.head = headkminusmer;
//                nodeL.end = endkminusmer;
//                nodeR.head = headkminusmer;
//                nodeR.end = endkminusmer;
//                nodeL.next.add(nodeR);
//                nodeR.prev.add(nodeL);

                // put into edge graph
                if (edges.containsKey(nodeL)) {
                    edges.get(nodeL).add(nodeR);
                } else {
                    HashSet<Node> set = new HashSet<Node>();
                    set.add(nodeR);
                    edges.put(nodeL, set);
                }
            }
        }

        printEdges(edges);
        HashSet<Node> heads = getHead(vertices, edges);

        // Find Contigs
        for (Node head : heads) {
            /*
            attach read to the end of contig
            for each read containing cur k-1:
                find cutoff
                substring of read
                if (contig contains substring of read)
                    curread = read
                    break
             */

            Node cur = head;
            String contig = "";
            String curRead = "";

            // Attach first read
            for (String read : cur.reads) {
                if (read.indexOf(cur.data) == 0) {
                    curRead = read;
                    contig += read;
                    //System.out.println(contig);
                    cur = vertices.get(read.substring(read.length() - k + 1));
                }
            }

            boolean endOfContig = false;
            boolean foundNext = false;

            while (!endOfContig) {
                foundNext = false;
                for (String read : cur.reads) {
                    if (!read.equals(curRead)){
                        // cutoff  is the index of the first character that
                        // comes immediately after the k-1 mer in read
                        int cutoff = read.indexOf(cur.data) + k - 1;
                        if (contig.contains(read.substring(0, cutoff))) {
                            foundNext = true;
                            curRead = read;
                            contig += read.substring(cutoff);
                            //System.out.println("read: "+read);
                            cur = vertices.get(read.substring(read.length() - k + 1));
                            //System.out.println(cur.data);
                            break;
                        }
                    }
                }
                if (!foundNext){
                    endOfContig = true;
                }
            }
            System.out.println(contig);
//            while (true){
////                System.out.println(contig);
////                System.out.println("prev: "+prev.data);
//                System.out.println("cur: "+cur.data);
//                if (edges.containsKey(cur)){
//                    for (int i = 0; i < cur.reads.size(); i ++){
//
//                        // requires k-1 < reads length
//                        if (cur.reads.get(i).indexOf(cur.data) == 0) {
//                            String str = cur.reads.get(i);
//                            contig += str;
//                            System.out.println(str +"hi");
//                            //prev = cur;
//                            cur = vertices.get(str.substring(str.length()-k+1, str.length()));
//                            break;
//                        } else {
//                            continue;
//                        }
//                    }
//                } else {
//                    break;
//                }
//            }
        }
    }


    private HashSet<Node> getHead(HashMap<String, Node> vertices, HashMap<Node, HashSet<Node>> edges){
        HashSet<Node> heads = new HashSet<Node>();
        for (Node vertex : vertices.values()) {
            boolean temp = true;
            for (HashSet<Node> listOfNode : edges.values()) {
                if (listOfNode.contains(vertex)){
                    temp = false;
                    continue;
                }

            }

            if (edges.containsKey(vertex) && temp){
                heads.add(vertex);
            }

        }
        System.out.println(heads.size()+"heads");

        for (Node node : heads) {
            System.out.println(node.data);
            System.out.println(node.reads.get(0));
        }
        System.out.println("\n");

        return heads;
    }



    /*assume a and b have equal length*/
    private double similarity(String a, String b) {
        if(a.length() == 0) return 1;
        int numberOfSimilarities = 0;
        for(int i = 0; i < a.length(); ++i) {
            if(a.charAt(i) == b.charAt(i)) {
                numberOfSimilarities++;
            }
        }
        return (double) numberOfSimilarities / a.length();
    }

    /*
     Help us visualize the edges of the de Bruijn edges
     */
    private void printEdges(HashMap<Node, HashSet<Node>> edges){
        for(Node nodeL:edges.keySet()) {
            for (Node node: edges.get(nodeL)){
                System.out.print(nodeL.data+" -> ");
                System.out.print(node.data+"\n");
//                System.out.print(nodeL.reads.get(0) + "(Lread)");
//                System.out.print(node.reads.get(0)+ "(Rread)");


            }
        }
        System.out.print("\n");
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
