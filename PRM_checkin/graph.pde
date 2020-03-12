//import java.util.Comparator; //<>//
//import java.util.Map;
//import java.util.Map.Entry;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;

class Graph{
  ArrayList<Vertex> vertexs=new ArrayList<Vertex>();
  ArrayList<Edge> edges=new ArrayList<Edge>();
  public void addVertex(Vertex vertex) {
    vertexs.add(vertex);
  }
  public void addEdge(Edge edge) {
      edges.add(edge);
      }
}

class Vertex{
  PVector pos_2D;
  boolean visited;
  Vertex current_pre;
  int index;
  Vertex(PVector l,int ind) {
    pos_2D = l.copy();
    index = ind;
    visited = false;
    }
  void draw_vert(){
    noStroke();
    fill(255, 0, 0);
    translate(pos_2D.x,cylinder_height,pos_2D.y);
    rotateX(PI/2);
    circle(0,0,5);
    rotateX(-PI/2);
    translate(-pos_2D.x,-cylinder_height,-pos_2D.y);
    fill(255);
  }
}

class Edge{
  Vertex start_pos;
  Vertex end_pos;
  float dist;
  Edge(Vertex start,Vertex end) {
    start_pos = start;
    end_pos = end;
    dist = PVector.dist(start_pos.pos_2D, end_pos.pos_2D);
    }
  void draw_edge(){
    stroke(5);
    fill(255, 0, 0);
    translate(0,cylinder_height,0);
    rotateX(PI/2);
    line(start_pos.pos_2D.x, start_pos.pos_2D.y, end_pos.pos_2D.x, end_pos.pos_2D.y);
    rotateX(-PI/2);
    translate(0,-cylinder_height,0);
    fill(255);
    noStroke();
  }
}


//void add_edge(){
//  for (int i = 0; i < vertexs.size(); i++) {
//    Vertex point1 = vertexs.get(i);
//    TreeMap<Integer,Float> distance_map = new TreeMap<Integer,Float>(new Comparator<Integer>() {
//      @Override
//      public int compare(Integer o1, Integer o2) {
//        return o2.compareTo(o1);
//      }
//    });
//    for (int j = 0; j < vertexs.size(); i++) {
//      if (i!=j){
//        Vertex point2 = vertexs.get(j);
//        float dist = point2.pos_2D.dist(point1.pos_2D); 
//        distance_map.put(j,dist);
//      }
//    }
//    Iterator<Entry<Integer, Float>> it = distance_map.entrySet().iterator();
//    for (int k = 0; k < 5; i++) {
//      Entry<Integer, Float> entry = it.next();
//      int closest_key = entry.getKey();
//      Edge new_edge =new Edge(point1,vertexs.get(closest_key));
//      graph.addEdge(new_edge);
//    }
//  }
//}

int[] indexesOfTopElements(float[] orig, int nummax) {
  float[] copy = Arrays.copyOf(orig,orig.length);
  Arrays.sort(copy);
  float[] honey = Arrays.copyOfRange(copy,copy.length - nummax, copy.length);
  int[] result = new int[nummax];
  int resultPos = 0;
  for(int i = 0; i < orig.length; i++) {
    float onTrial = orig[i];
    int index = Arrays.binarySearch(honey,onTrial);
    if(index < 0) continue;
    result[resultPos++] = i;
  }
  return result;
}

ArrayList<Vertex> BFS(Graph graph) {
  ArrayList<Vertex> solution = new ArrayList<Vertex>();
  ArrayList<Vertex> vertexs=graph.vertexs;
  ArrayList<Edge> edges=graph.edges;
  Queue<Vertex> queue = new LinkedList<Vertex>();
  //add start point to the queue
  queue.add(vertexs.get(0));
  vertexs.get(0).visited=true;
  //System.out.print(vertexs.get(0));
  
  while(!queue.isEmpty()) {
    Vertex vertex=queue.peek();
    //println(vertex.index);
    //println(vertex.pos_2D);
    //println(queue.size());
    for(Edge edge:edges) {
      if(edge.start_pos.equals(vertex)&&edge.end_pos.visited==false) {
        queue.add(edge.end_pos);
        edge.end_pos.visited=true;
        edge.end_pos.current_pre = edge.start_pos;
        //System.out.print(edge.end_pos);
        if (edge.end_pos == vertexs.get(1)){
          Vertex vertex0 = vertex;
          println("find!");
          while(vertex0.current_pre != null){
            solution.add(vertex0.current_pre);
            //print(vertex0.pos_2D);
            vertex0 = vertex0.current_pre;
          }
          return solution;
        }
      }
    }
    queue.poll();
  }
  return solution;
}
