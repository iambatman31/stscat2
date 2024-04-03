//BOUNDARY TRAVERSAL OF BINARY TREE

class TreeNode {
    int val;
    TreeNode left, right;

    public TreeNode(int value) {
        val = value;
        left = right = null;
    }
}

class BoundaryTraversal {

    private static void printLeaves(TreeNode node) {
        if (node != null) {
            printLeaves(node.left);
            if (node.left == null && node.right == null) System.out.print(node.val + " ");
            printLeaves(node.right);
        }
    }

    private static void printLeftBoundary(TreeNode node) {
        if (node != null) {
            if (node.left != null) {
                System.out.print(node.val + " ");
                printLeftBoundary(node.left);
            } else if (node.right != null) {
                System.out.print(node.val + " ");
                printLeftBoundary(node.right);
            }
        }
    }

    private static void printRightBoundary(TreeNode node) {
        if (node != null) {
            if (node.right != null) {
                printRightBoundary(node.right);
                System.out.print(node.val + " ");
            } else if (node.left != null) {
                printRightBoundary(node.left);
                System.out.print(node.val + " ");
            }
        }
    }

    public static void boundaryTraversal(TreeNode root) {
        if (root != null) {
            System.out.print(root.val + " "); 
            printLeftBoundary(root.left);

            printLeaves(root.left);
            printLeaves(root.right);

            printRightBoundary(root.right);
        }
    }

    public static void main(String[] args) {
        TreeNode root = new TreeNode(8);
        root.left = new TreeNode(3);
        root.right = new TreeNode(10);
        root.left.left = new TreeNode(1);
        root.left.right = new TreeNode(6);
        root.left.right.left = new TreeNode(4);
        root.left.right.right = new TreeNode(7);
        root.right.right = new TreeNode(14);
        root.right.right.left = new TreeNode(13);

        boundaryTraversal(root);
    }
}

//BFS

import java.util.*;

class Graph {
  private int V; 
  private LinkedList<Integer>[] adj; 


  Graph(int v) {
    V = v;
    adj = new LinkedList[v];
    for (int i = 0; i < v; ++i) {
      adj[i] = new LinkedList<>();
    }
  }


  void addEdge(int v, int w) {
    adj[v].add(w);
  }


  void BFS(int s) {
    boolean[] visited = new boolean[V];
    Queue<Integer> queue = new LinkedList<>();
    visited[s] = true;
    queue.add(s);

    while (!queue.isEmpty()) {
      s = queue.poll();
      System.out.print(s + " ");
      Iterator<Integer> i = adj[s].listIterator();
      while (i.hasNext()) {
        int n = i.next();
        if (!visited[n]) {
          visited[n] = true;
          queue.add(n);
        }
      }
    }
  }

  public static void main(String[] args) {
    Graph g = new Graph(4);

    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 0);
    g.addEdge(2, 3);
    g.addEdge(3, 3);

    System.out.println("Following is Breadth First Traversal (starting from vertex 2)");
    g.BFS(2);
  }
}

//DFS

import java.util.*;
class Graph {
  private LinkedList<Integer> adjLists[];
  private boolean visited[];
  Graph(int vertices) {
    adjLists = new LinkedList[vertices];
    visited = new boolean[vertices];
    for (int i = 0; i < vertices; i++)
      adjLists[i] = new LinkedList<Integer>();
  }
void addEdge(int src, int dest) {
    adjLists[src].add(dest);
  }// DFS algorithm
  void DFS(int vertex) {
    visited[vertex] = true;
System.out.print(vertex + " ");
    Iterator<Integer> ite = adjLists[vertex].listIterator();
    while (ite.hasNext()) {
      int adj = ite.next();
      if (!visited[adj])
        DFS(adj);
    }  }
public static void main(String args[]) {
    Graph g = new Graph(4);

    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);

    System.out.println("Following is Depth First Traversal");

    g.DFS(2);
  }
}

//BST

public class BST {
    TreeNode firstIncorrectNode = null;
    TreeNode secondIncorrectNode = null;
    TreeNode prevNode = new TreeNode(Integer.MIN_VALUE);

    public void recoverTree(TreeNode root) {
        inorder(root);
        int temp = firstIncorrectNode.val;
        firstIncorrectNode.val = secondIncorrectNode.val;
        secondIncorrectNode.val = temp;
    }

    private void inorder(TreeNode node) {
        if (node == null) return;
        inorder(node.left);
        if (firstIncorrectNode == null && prevNode.val >= node.val) {
            firstIncorrectNode = prevNode;
        }
        if (firstIncorrectNode != null && prevNode.val >= node.val) {
            secondIncorrectNode = node;
        }
        prevNode = node;
        inorder(node.right);
    }

    public static void main(String[] args) {
        TreeNode root = new TreeNode(3);
        root.left = new TreeNode(1);
        root.right = new TreeNode(4);
        root.right.left = new TreeNode(2);
        Main solution = new Main();
        solution.recoverTree(root);
        System.out.println("Inorder Traversal of Recovered BST:");
        printInorder(root);
    }

    private static void printInorder(TreeNode node) {
        if (node == null) return;
        printInorder(node.left);
        System.out.print(node.val + " ");
        printInorder(node.right);
    }
}

//VERTICAL ORDER TRAVERSAL

import java.util.*;
import java.util.AbstractMap.SimpleEntry;
class TreeNode {
    int val;
    TreeNode left;
    TreeNode right;
    public TreeNode(int val) {
        this.val = val;
left = null;
        right = null;
    }
}
public class vertical_order_traversal {

    public static List<List<Integer>>
verticalOrderTraversal(TreeNode root) {
        List<List<Integer>> result = new
ArrayList<>();
if (root == null) {
            return result;
}
        Map<Integer, List<Integer>>
verticalMap = new TreeMap<>();
        Queue<SimpleEntry<TreeNode,
Integer>> nodeQueue = new LinkedList<>();
        nodeQueue.offer(new
SimpleEntry<>(root, 0));
        while (!nodeQueue.isEmpty()) {
            SimpleEntry<TreeNode,
Integer> entry = nodeQueue.poll();
            TreeNode node =
entry.getKey();
            int col = entry.getValue();
verticalMap.computeIfAbsent(col, k -> new
ArrayList<>()).add(node.val);
            if (node.left != null) {
                nodeQueue.offer(new
SimpleEntry<>(node.left, col - 1));
}
if (node.right != null) {
    nodeQueue.offer(new
SimpleEntry<>(node.right, col + 1));
}
}
for (List<Integer> values :
verticalMap.values()) {
result.add(values);
}
return result;
}
public static void main(String[] args) {
        // Sample binary tree input
        TreeNode root = new TreeNode(1);
        root.left = new TreeNode(2);
        root.right = new TreeNode(3);
        root.left.left = new TreeNode(4);
        root.left.right = new TreeNode(5);
        root.right.left = new TreeNode(6);
        root.right.right = new TreeNode(7);
        List<List<Integer>>
        verticalOrderResult =
        verticalOrderTraversal(root);

        for (List<Integer> column :
        verticalOrderResult) {
        for (int val : column) {
        System.out.print(val + " ");
        }
        System.out.println();
}
} 
}
return horizontalView;
    }
    public static void main(String[]args) {
        // Sample binary tree input
        TreeNode root = new
TreeNode('A');
        root.left = new TreeNode('B');
        root.right = new TreeNode('C');
        root.left.left = new
TreeNode('D');
        root.left.right = new
TreeNode('E');
        root.right.left = new
TreeNode('F');
        root.right.right = new
TreeNode('G');
List<Character> horizontalViewResult =
horizontalView(root);
        // Printing the Horizontal View
        System.out.print("HorizontalView: ");
        for (char node :
horizontalViewResult) {
            System.out.print(node + " ");
        }
        System.out.println();
    }
}

//VIEW OF TREE

import java.util.*;
class TreeNode {
    char val;
    TreeNode left;
    TreeNode right;
    public TreeNode(char val) {
        this.val = val;
left = null;
        right = null;
    }}
public class tree {
    public static List<Character>
horizontalView(TreeNode root) {
        List<Character> horizontalView =
new ArrayList<>();
        if (root == null) {
            return horizontalView;
        }
Queue<TreeNode> queue = new
LinkedList<>();
        queue.offer(root);
        while (!queue.isEmpty()) {
            int levelSize = queue.size();
            for (int i = 0; i <
levelSize; i++) {
                TreeNode node =
queue.poll();
            horizontalView.add(node.val);
                if (node.left != null) {
queue.offer(node.left);
                }
                if (node.right != null) {
queue.offer(node.right);
                }
} }
return horizontalView;
    }
    public static void main(String[]
args) {
        // Sample binary tree input
        TreeNode root = new
TreeNode('A');
        root.left = new TreeNode('B');
        root.right = new TreeNode('C');
        root.left.left = new
TreeNode('D');
        root.left.right = new
TreeNode('E');
        root.right.left = new
TreeNode('F');
        root.right.right = new
TreeNode('G');
List<Character> horizontalViewResult =
horizontalView(root);
        // Printing the Horizontal View
        System.out.print("Horizontal View: ");
        for (char node :horizontalViewResult) {
            System.out.print(node + " ");
        }
        System.out.println();
    }
}

//DIAL'S ALGORITHM

import java.util.*;
class Graph {
    private int V;
    private List<List<Node>> adj;

    public Graph(int V) {
        this.V = V;
        adj = new ArrayList<>(V);
        for (int i = 0; i < V; i++) {
            adj.add(new ArrayList<>());
        }
    }
    public void addEdge(int source, int destination, int weight) {
        Node node = new Node(destination, weight);
        adj.get(source).add(node);
    }
    public void dijkstra(int startVertex) {
 int[] distance = new int[V];
        Arrays.fill(distance, Integer.MAX_VALUE);

        distance[startVertex] = 0;

        PriorityQueue<Node> pq = new PriorityQueue<>(V, Comparator.comparingInt(node -> node.weight));
        pq.add(new Node(startVertex, 0));
while (!pq.isEmpty()) {
            int currentVertex = pq.poll().vertex;
            for (Node neighbor : adj.get(currentVertex)) {
                int newDist = distance[currentVertex] + neighbor.weight;
               
if (newDist < distance[neighbor.vertex]) {
          distance[neighbor.vertex] = newDist;
             pq.add(new Node(neighbor.vertex, newDist));
                        }
                    }
                }
        // Print the distances
        System.out.println("Vertex\tDistance from Source");
                for (int i = 0; i < V; i++) {
                    System.out.println(i + "\t" + distance[i]);
                }
            }
        static class Node {
                int vertex;
                int weight;
             public Node(int vertex, int weight) {  
        this.vertex = vertex;
                    this.weight = weight;
                }
            }
        }
        public class Main {
            public static void main(String[] args) {
                int V = 5; // Number of vertices
                int source = 0; // Source vertex
                Graph graph = new Graph(V);
                graph.addEdge(0, 1, 2);
                graph.addEdge(0, 3, 1);
                graph.addEdge(1, 2, 3);
                graph.addEdge(1, 3, 2);
                graph.addEdge(3, 4, 4);
                graph.addEdge(4, 2, 1);
        
                graph.dijkstra(source);
            }
        }

//BELLMAN FORD ALGORITHM

public class Main {
        static class CreateGraph {
            class CreateEdge {
                int src, dest, weight;
                CreateEdge(int src, int dest, int weight) {
                    this.src = src;
                    this.dest = dest;
                    this.weight = weight;
                }
            }
            int V, E;
            CreateEdge edge[];
            CreateGraph(int v, int e) {
                V = v;
                E = e;
                edge = new CreateEdge[e];
            }
           
    void BellmanFord(int src) {
                int dist[] = new int[V];
                for (int i = 0; i < V; ++i)
                dist[i] = Integer.MAX_VALUE;
                dist[src] = 0;
       for (int i = 1; i < V; ++i) {
        for (int j = 0; j < E; ++j) {
          int u = edge[j].src;
             int v = edge[j].dest;
                int w = edge[j].weight;
                 if (dist[u] != Integer.MAX_VALUE && dist[u] + w < dist[v])
                dist[v] = dist[u] + w;
                    }
                }
          for (int j = 0; j < E; ++j) {
                    int u = edge[j].src;
                    int v = edge[j].dest;
    int w = edge[j].weight;
              if (dist[u] != Integer.MAX_VALUE && dist[u] + w < dist[v]) {                   System.out.println("Graph contains negative weight cycle");
                    return;
                    }
                }
                printSolution(dist);
            }
            void printSolution(int dist[]) {
                System.out.println("Vertex Distance from Source");
                for (int i = 0; i < V; ++i)
    System.out.println(i + "\t\t" + dist[i]);
            }
        }
    public static void main(String[] args) {
            int V = 5;
            int E = 7;
    CreateGraph graph = new CreateGraph(V, E);
            graph.edge[0] = graph.new CreateEdge(0, 1, 5);
            graph.edge[1] = graph.new CreateEdge(0, 2, 4);
            graph.edge[2] = graph.new CreateEdge(1, 3, 3);
            graph.edge[3] = graph.new CreateEdge(2, 1, 6);
            graph.edge[4] = graph.new CreateEdge(3, 2, 2);
            graph.edge[5] = graph.new CreateEdge(1, 4, -4);
            graph.edge[6] = graph.new CreateEdge(4, 2, 2);
            graph.BellmanFord(0); // 0 is the source vertex
        }}


//TOPOLOGICAL SORT

import java.util.*;
class Graph {
    private Map<Integer, List<Integer>> adjacencyList;
    private int vertices;
    public Graph(int vertices) {
        this.vertices = vertices;
        this.adjacencyList = new HashMap<>();
        for (int i = 0; i < vertices; i++) {
            this.adjacencyList.put(i, new ArrayList<>());
        }
    }
    public void createEdge(int u, int v) {
        this.adjacencyList.get(u).add(v);
    }
    public void topologicalSort() {
        int[] totalIndegree = new int[vertices];
       for (int i = 0; i < vertices; i++) {
            for (int j : adjacencyList.get(i)) {
                totalIndegree[j]++;
            }}
Queue<Integer> queue = new LinkedList<>();
        for (int i = 0; i < vertices; i++) {
            if (totalIndegree[i] == 0) {
                queue.add(i);
            }
        }
        int visitedNodes = 0;
        List<Integer> order = new ArrayList<>();

        while (!queue.isEmpty()) {
            int u = queue.poll();
            order.add(u);

            for (int i : adjacencyList.get(u)) {
                totalIndegree[i]--;

                if (totalIndegree[i] == 0) {
                    queue.add(i);
                }
            }
            visitedNodes++;
        }
if (visitedNodes != vertices) {
    System.out.println("There's a cycle present in the Graph.\nGiven graph is not DAG");
} else {
    System.out.println(order);
}
}

public static void main(String[] args) {
Graph graph = new Graph(6);
graph.createEdge(0, 1);
graph.createEdge(0, 2);
graph.createEdge(1, 3);
graph.createEdge(1, 5);
graph.createEdge(2, 3);
graph.createEdge(2, 5);
graph.createEdge(3, 4);
graph.createEdge(5, 4);

graph.topologicalSort();
}
}

// HEAP SORT

class HeapSort {
        static void heapify(int a[], int n, int i) {
            int largest = i;
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            if (left < n && a[left] > a[largest])
                largest = left;
            if (right < n && a[right] > a[largest])
                largest = right;
            if (largest != i) {
                int temp = a[i];
                a[i] = a[largest];
                a[largest] = temp;
                heapify(a, n, largest);
            }
        }
     static void heapSort(int a[], int n) {
      for (int i = n / 2 - 1; i >= 0; i--)
                heapify(a, n, i);
    
            for (int i = n - 1; i >= 0; i--) {
                int temp = a[0];
                a[0] = a[i];
                a[i] = temp;
    
                heapify(a, i, 0);
            }
        }
    
        static void printArr(int a[], int n) {
    
        for (int i = 0; i < n; ++i)
        System.out.print(a[i] + " ");
}

public static void main(String args[]) {
    int a[] = {35, 17, 10, 90, 24, -3, -8};
    int n = a.length;

    System.out.println("Original Array:");
    printArr(a, n);

    heapSort(a, n);

    System.out.println("\nSorted Array:");
    printArr(a, n);
}
}


// BINOMIAL HEAP

import java.io.*;
class BinomialHeapNode {
  int key, degree;
  BinomialHeapNode parent;
  BinomialHeapNode sibling;
  BinomialHeapNode child;
  public BinomialHeapNode(int k)
  {

  key = k;
  degree = 0;
  parent = null;
  sibling = null;
  child = null;
  }
  public BinomialHeapNode reverse(BinomialHeapNode sibl)
  {
 
BinomialHeapNode ret;
  if (sibling != null)
  ret = sibling.reverse(this);
  else
  ret = this;
  sibling = sibl;
  return ret;
  }
  public BinomialHeapNode findMinNode()
  {
  BinomialHeapNode x = this, y = this;
  int min = x.key;
  while (x != null) {
 
    if (x.key < min) {
          y = x;
          min = x.key;
          }
        
          x = x.sibling;
          }
          return y;
          }
        public BinomialHeapNode findANodeWithKey(int value)
          {
          BinomialHeapNode temp = this, node = null;
          while (temp != null) {
          if (temp.key == value) {
          node = temp;
          break;
          }
          if (temp.child == null)
          temp = temp.sibling;
          else {
          node=temp.child.findANodeWithKey(value);
          if (node == null)
             temp = temp.sibling;
          else
          break;
          }
          }
          return node;
          }
          public int getSize()
          {
          return (1 + ((child == null) ? 0 : child.getSize())+ ((sibling == null) ? 0 : sibling.getSize()));
          }
        }
        class BinomialHeap {
          private BinomialHeapNode Nodes;
          private int size;
          public BinomialHeap()
          {
          Nodes = null;
          size = 0;
          }
          public boolean isEmpty() { return Nodes == null; }
          public int getSize() { return size; }
          public void makeEmpty()
          {
          Nodes = null;
          size = 0;
          }
          public void insert(int value)
          {
        
          if (value > 0) {
        BinomialHeapNode temp= new BinomialHeapNode(value);
          if (Nodes == null) {
                Nodes = temp;
           size = 1;
          }
        else {
          unionNodes(temp);
          size++;
          }
          }
          }
          private void merge(BinomialHeapNode binHeap)
          {
          BinomialHeapNode temp1 = Nodes, temp2 = binHeap;
          while ((temp1 != null) && (temp2 != null)) {
          if (temp1.degree == temp2.degree) {
         
            BinomialHeapNode tmp = temp2;
              temp2 = temp2.sibling;
              tmp.sibling = temp1.sibling;  temp1.sibling = tmp;
              temp1 = tmp.sibling;
              }
            else {
              if (temp1.degree < temp2.degree) {
              if ((temp1.sibling == null)|| (temp1.sibling.degree> temp2.degree)) {  BinomialHeapNode tmp = temp2;  temp2 = temp2.sibling;
              tmp.sibling = temp1.sibling;
              temp1.sibling = tmp;
              temp1 = tmp.sibling;
              }
              else {  temp1 = temp1.sibling;
              }
              }
            else {
            BinomialHeapNode tmp = temp1;
              temp1 = temp2;
              temp2 = temp2.sibling;
              temp1.sibling = tmp;
              if (tmp == Nodes) {
              Nodes = temp1;
              }
              else {
              }
              }
              }
              }
              if (temp1 == null) {
              temp1 = Nodes;
              while (temp1.sibling != null) {
              temp1 = temp1.sibling;
              }
             
public class main {
    public static void main(String[] args)
        {
            BinomialHeap binHeap = new BinomialHeap();
            binHeap.insert(12);
  binHeap.insert(8);
  binHeap.insert(5);
  binHeap.insert(15);
  binHeap.insert(7);
  binHeap.insert(2);
  binHeap.insert(9);
  System.out.println("Size of the binomial heap is "+ binHeap.getSize());
  binHeap.displayHeap();
 
binHeap.delete(15);
  binHeap.delete(8);
  System.out.println("Size of the binomial heap is "+ binHeap.getSize());
  binHeap.displayHeap();
  binHeap.makeEmpty();
  System.out.println(binHeap.isEmpty());
  }
}


// WINNER TREE

public class WinnerTree {
        private int[] tree;
        private int[] players;
        private int numPlayers;
    
        public WinnerTree(int[] players) {
            this.players = players;
            this.numPlayers = players.length;
            int treeSize = 2 * numPlayers - 1;
            this.tree = new int[treeSize];
            buildWinnerTree();
        }
    
        private void buildWinnerTree() {
            for (int i = numPlayers - 1; i < tree.length; i++) {
                tree[i] = i - (numPlayers - 1);
            }
           
    for (int i = numPlayers - 2; i >= 0; i--) {
                tree[i] = players[tree[2 * i + 1]] < players[tree[2 * i + 2]] ? tree[2 * i + 2] : tree[2 * i + 1];
            }
        }
    public int getWinner() {
            return players[tree[0]];
        }
    
        public static void main(String[] args) {
            int[] players = {3,7,1,9,4,2,8,5};
            WinnerTree winnerTree = new WinnerTree(players);
            System.out.println("The winner is: " + winnerTree.getWinner());
        }
    }
    
    
// WINNER TREE SCORE

import java.util.Arrays;
public class MinWinnerTree {
    private int[] tree;
    private int[] players;
public MinWinnerTree(int[] players) {
        this.players = players;
        int n = players.length;
        int treeSize = calculateTreeSize(n);
        tree = new int[2 * treeSize - 1];
        Arrays.fill(tree, -1);
        for (int i = 0; i < n; i++) {
            tree[treeSize - 1 + i] = i;
        }
        buildWinnerTree(0, 0, treeSize - 1);
    }
private int calculateTreeSize(int n) {
        int treeSize = 1;
        while (treeSize < n) {
            treeSize *= 2;
        }
        return treeSize;
    }
private void buildWinnerTree(int node, int left, int right) {
        if (left == right) return;
        int mid = (left + right) / 2;
        buildWinnerTree(2 * node + 1, left, mid);
        buildWinnerTree(2 * node + 2, mid + 1, right);
        tree[node] = players[tree[2 * node + 1]] <= players[tree[2 * node + 2]] ? tree[2 * node + 1] : tree[2 * node + 2];
    }
public int getWinnerIndex() {
        return tree[0];
    }
public static void main(String[] args) {
        int[] players = {3, 7, 1, 9, 4, 2, 8, 5};
        MinWinnerTree winnerTree = new MinWinnerTree(players);
        int winnerIndex = winnerTree.getWinnerIndex();
        int winningScore = players[winnerIndex];
        System.out.println("The player with the lowest score is at index: " + winnerIndex);
        System.out.println("The score of the winning player is: " + winningScore);
    }
}

                                                    

