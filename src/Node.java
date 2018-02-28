import java.util.ArrayList;
import java.util.HashSet;

public class Node {
    public String data;
    public ArrayList<String> reads = new ArrayList<String>();
    public String head;
    public String end;
    public boolean isVisited = false;
    public ArrayList<Node> prev = new ArrayList<>();
    public ArrayList<Node> next = new ArrayList<>();

    public Node(String data, String read){

        this.data = data;
        if (!reads.contains(read)){
            reads.add(read);
        }


    }

    public void addReads(String read){
        if (!reads.contains(read)){
            reads.add(read);
        }
    }
    // might implement our own hashcode function
}
