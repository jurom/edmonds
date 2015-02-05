import 'dart:io';

/*
 * Warning: Running this code may output debug prints to file 'out.txt' if debug == true.
 * 
 * General structure:
 * 
 * There exists an instance of Graph, which contains all the necessary information about the given graph (edges, weights),
 * converts the input into some different objects (Flower-s, Edge-s.. ) and eventually runs and resolves all the iterations.
 * Additionally, all pairing (blocking) edges are all time in a global variable edgesM (edgesL).
 * 
 * As for the data structures, there exists a Flower, which can be either a ComplexFlower (containing other flowers) or SingleVertex 
 * (single, blue flower). As those two types of flowers had slightly different behavior, I found it more convenient to make them
 * inherit some common implementation from abstract Flower and then implement something themselves. These Flower-s are able to
 * determine the max charge that may be added (subtracted) by themselves, as well as recursively find an alternating path starting in
 * a given vertex in peduncle and ending in arbitrary (almost) vertex. A flower remembers its wrapping flower, and actually the
 * whole list of wrapping flowers (this was formerly always recalculated by getters, but for slight optimisation for larger inputs
 * it is not recalculated only when necessary). Aside from its wrapping flower, the outermost flower remembers the parentStructure it 
 * is in (Dumbbell or HungarianTreeNode).
 * Another data structure is an Edge. An Edge is capable of determining the maximum charge, which may be added to it. It remembers the
 * vertices it connects, an therefore has full information about their wrapping flowers and parent structures - therefore the edge alone
 * is able to determine, what are the two structures it connects and find the proper max charge.
 * Then, there are data structures wrapping the Flower-s. One of them is a Dumbbell - it's truly dumb, it basically only remembers the
 * flowers and pairing edge (which I have acutally never used). And the other one is a HungarianTreeNode (HTN). This HTN remembers what
 * a regular node in double referenced tree should remember - children and parent. HTN knows its depthParity (0 or 1) and the max charge
 * to add (subtract). Of course, as there exists a node, we must have a tree as well - HungarianTree. It is basically just capable of 
 * finding the max charge and adding charge recursively.
 * 
 * Workflow:
 * 
 * An instance of Graph constructs SingleVertex from every given vertex, and Edge for every given edge. All the time, graph keeps
 * references to all the trees and dumbbells. Every iteration, the graph finds out the maximum charge able to add to every tree, 
 * adds the charge, and resolves the aroused issues (using standard algorithm for P1,.. P4). It tries to resolve all of them, even though
 * some of them may not be issues anymore. 
 * Additionally, if an edge emits the maximum charge able to add larger than maximum weight, this means that any charge may be 
 * added - therefore, if the maximal charge amongst all the trees is larger than maximum and pairing edges do not form a perfect pairing yet,
 * there is no perfect matchin in this graph.
 *  
 */



/// Pairing edges
Set<Edge> edgesM = new Set();
/// Blocking edges
Set<Edge> edgesL = new Set();
num maxWeight = 0;
infiniteCharge() => maxWeight+1;

// Just for calculation time measuring purposes
Stopwatch epsTime = new Stopwatch();
List<List> solutionTime = new Iterable.generate(4,(i) => [new Stopwatch(), 0]).toList();
start(int i) {
  solutionTime[i][0].start();
  solutionTime[i][1]++;
}
stop(i) => solutionTime[i][0].stop();

// All potential problems are every iteration stored here under key equal to the added (subtracted) charge that would cause them
Map<num, Queue> potentialProblems = {};

/// Add a potential problem when checking for max charge
addPotentialProblem(dynamic problem, num eps) {
  if (potentialProblems.containsKey(eps)) {
    potentialProblems[eps].put(problem);
  } else {
    potentialProblems[eps] = new Queue()..put(problem);
  }
  return eps;
}

getProblems(num eps) => potentialProblems[eps];

clearProblems() => potentialProblems.clear();

class Queue {

  Set _q = new Set();

  Queue();

  put(dynamic a) => _q.add(a);
  get() {
    dynamic a = _q.first;
    _q.remove(a);
    return a;
  }

  get isEmpty => _q.isEmpty;
  get isNotEmpty => _q.isNotEmpty;

  toString() => _q.toString();

}

min(Iterable<Comparable> l) {
  if (l.length == 0) return null;
  var m = l.first;
  l.forEach((e) => e < m ? m = e : null);
  return m;
}

concat(List<List> lists) {
  List conc = [];
  lists.forEach((List l) => conc.addAll(l));
  return conc;
}

/// Blocking edges in edgeborder of [f]
getBlockingEdges(Flower f) => f.edgeBorder.intersection(edgesL);
/// Blocking edges in edgeborder of [f], which do not leave the wrapping flower (useful in paths)
getInnerBlockingEdges(Flower f) => f.innerEdgeBorder.intersection(edgesL);

/// Returns flowers on a desired (same) level for given [edge].
List<Flower> getEdgeFlowers(Edge edge, num levelFromTop) {
  List l1 = edge.v1.wrappingFlowerList;
  List l2 = edge.v2.wrappingFlowerList;
  List edgeFlowers = [l1[l1.length - 1 - levelFromTop], l2[l2.length - 1 - levelFromTop]];
  return edgeFlowers;
}

/// Finds an alternating path of flowers on the same level of even length from peduncle to a flower that satisfied given function
List<Flower> findEvenPath(Flower peduncle, bool verifyDesiredEnd(Flower f), num levelFromTop) {

  if (verifyDesiredEnd(peduncle)) return [peduncle];

  Flower getOtherFlower(Edge edge, Flower flower) => getEdgeFlowers(edge, levelFromTop).where((Flower f) => !identical(f, flower)).single;

  for (Edge edge in getInnerBlockingEdges(peduncle)) {
    Flower current = peduncle;
    List path = [current];
    current = getOtherFlower(edge, current);
    while (!verifyDesiredEnd(current)) {
      path.add(current);
      // Use edge from M (L resp.), which is the outermost inner edge for the enclosing flower
      current = getOtherFlower(current.innerEdgeBorder.intersection(path.length % 2 == 0 ? edgesM : edgesL).single, current);
    }
    if (path.length % 2 == 0) {
      // Found the proper path (odd vertices, even length)
      path.add(current);
      return path;
    }
  }


  // As there is always odd number of Flower-s in a ComplexFlower, this point is never reached
  assert(false);
  // Only for syntax highliting..
  return [];
}


class Edge {

  num weight;
  SingleVertex v1;
  SingleVertex v2;

  List<SingleVertex> get vertices => [v1, v2];

  Set<Flower> get crossedBy => v1.wrappingFlowers.union(v2.wrappingFlowers);

  get charge => crossedBy.map((Flower f) => f.charge).reduce((x, y) => x + y);

  bool get checkSkip => ((charge != weight) || edgesM.contains(this) || edgesL.contains(this));
  bool get checkP2 {
    if (checkSkip) return false;
    return (((v1.parentStructure is Dumbbell) && (v2.parentStructure is HungarianTreeNode) && (v2.parentStructure.depthParity == 0)) || (v2.parentStructure is Dumbbell) && (v1.parentStructure is HungarianTreeNode) && (v1.parentStructure.depthParity == 0));
  }

  bool get checkP3 {
    if (checkSkip) return false;
    return ((v1.parentStructure is HungarianTreeNode) && (v2.parentStructure is HungarianTreeNode)
        && (identical(v1.parentStructure.root, v2.parentStructure.root))
        && (!identical(v1.wrappingFlowerList.last, v2.wrappingFlowerList.last))
        && (v1.parentStructure.depthParity + v2.parentStructure.depthParity == 0));
  }

  bool get checkP4 {
    if (checkSkip) return false;
    return ((v1.parentStructure is HungarianTreeNode) && (v2.parentStructure is HungarianTreeNode) && (!identical(v1.parentStructure.root, v2.parentStructure.root)) && (v1.parentStructure.depthParity + v2.parentStructure.depthParity == 0));
  }

  num get maxChargeToAdd {
    if (edgesL.contains(this) || edgesM.contains(this)) return infiniteCharge();
    dynamic parentStruct1 = v1.parentStructure;
    dynamic parentStruct2 = v2.parentStructure;

    // In a tree, on levels with different parity
    if ((parentStruct1 is HungarianTreeNode) && (parentStruct2 is HungarianTreeNode) && (parentStruct1.depthParity + parentStruct2.depthParity == 1)) {
//      printF("Double tree, unequal parity: ${infiniteCharge()}");
      return infiniteCharge();
    }

    // In a tree, on levels with equal (even) parity
    if (((parentStruct1 is HungarianTreeNode) && (parentStruct2 is HungarianTreeNode) && (parentStruct1.depthParity + parentStruct2.depthParity == 0))) {
//      printF("Double tree, equal parity: $weight, $charge, ${(weight - charge)/2}");
      return addPotentialProblem(this, (weight - charge) / 2);
    }

    // Tree + dumbbell
    if ((parentStruct1 is HungarianTreeNode) && (parentStruct2 is Dumbbell) && (parentStruct1.depthParity == 0)) {
//      printF("Tree + dumbbell: ${(weight - charge)}");
      return addPotentialProblem(this, (weight - charge));
    }

    // Dumbbell + tree
    if ((parentStruct2 is HungarianTreeNode) && (parentStruct1 is Dumbbell) && (parentStruct2.depthParity == 0)) {
//      printF("Tree + dumbbell: ${(weight - charge)}");
      return addPotentialProblem(this, (weight - charge));
    }

    // Otherwise, possibilities are endless..
    return infiniteCharge();
  }

  Edge(SingleVertex _v1, SingleVertex _v2, this.weight) {
    if (_v1.vertex < _v2.vertex) {
      v1 = _v1;
      v2 = _v2;
    } else {
      v1 = _v2;
      v2 = _v1;
    }
  }

  toString() => "${v1.vertex} ${v2.vertex} C: $charge";

}

class Graph {
  Map _weights;
  num verticesCount;
  Set<HungarianTree> trees = new Set();
  Set<Dumbbell> dumbbells = new Set();
  Map<num, SingleVertex> numberToVertex = {};

  static Map<num, Set<Edge>> edges = {};

  static Set getNeighbors(num a) => edges[a];

  Graph(this.verticesCount, List<List> _edges) {
    addEdge(num v, Edge e) {
      if (edges.containsKey(v)) {
        edges[v].add(e);
      } else {
        edges[v] = new Set.from([e]);
      }
    }
    _edges.forEach((e) {
      if (e[2] > maxWeight) maxWeight = e[2];
      [e[0], e[1]].forEach((num e) {
        if (!numberToVertex.containsKey(e)) {
          numberToVertex[e] = new SingleVertex(e);
          HungarianTreeNode htn = new HungarianTreeNode(numberToVertex[e], null, 0);
          numberToVertex[e].parentStructure = htn;
          HungarianTree ht = new HungarianTree(htn);
          htn.tree = ht;
          trees.add(ht);
        }
      });
      Edge ed = new Edge(numberToVertex[e[0]], numberToVertex[e[1]], e[2]);
      addEdge(e[0], ed);
      addEdge(e[1], ed);
    });
  }

  // Assert - every vertex is always incident with at most one edge from M
  get isPerfectPairing => (verticesCount % 2 == 0) && (edgesM.length == verticesCount ~/ 2);
  
  bool _no1Factor = false;
  
  get no1Factor => _no1Factor;

  nextIteration() {
    if (isPerfectPairing || no1Factor) return;
    Stopwatch tmpEps = new Stopwatch();
    tmpEps.start();
    epsTime.start();
    clearProblems();
    List charges = trees.map((HungarianTree ht) => ht.maxChargeToAdd).toList();
    num maxChargeToAdd = min(charges);
    printF("Max charge to add: $maxChargeToAdd");
    
    if (maxChargeToAdd == infiniteCharge()) {
      // No problem will be aroused by adding any charge - it does not have 1 factor
      _no1Factor = true;
      return;
    }
    trees.forEach((HungarianTree ht) => ht.addCharge(maxChargeToAdd));
    tmpEps.stop();
    epsTime.stop();

    Queue problems = getProblems(maxChargeToAdd);

    print("Found eps in ${tmpEps.elapsed}");
    String str = "";
    numberToVertex.forEach((k, v) => str += "$k:${v.charge}, ");
    while(problems.isNotEmpty) {
      printF("Pairing: $edgesM");
      printF("Blocking: $edgesL");
      printF("Number of trees: ${trees.length}");
      printF("Trees: \n ${trees}");
      printF("Dumbbells: \n ${dumbbells}");

      dynamic problem = problems.get();
      printF("Current problem: $problem");
      if (problem is ComplexFlower) {
        start(0);
        // P1
        printF("P1");
        // It may not be a problem anymore (fixed by something else)
        if ((problem.charge != 0) || (problem.wrappingFlower != null) || (problem.parentStructure is! HungarianTreeNode)
            || (problem.parentStructure.depthParity == 0)) return;
        assert(problem.parentStructure is HungarianTreeNode); // for code completion suggestions

        // Debug print
        printF("Edge border of problem: ${problem.edgeBorder.map((Edge e) {
          if (edgesL.contains(e)) {
            assert(!edgesM.contains(e));
            return "${e}L";
          } else if (edgesM.contains(e)) {
            return "${e}M";
          }
          return e.toString();
        }).toList()}");

        Edge blockingEdge = problem.edgeBorder.intersection(edgesL).single;

        SingleVertex topFlowerVertex = blockingEdge.vertices.where((SingleVertex v) => v.wrappingFlowers.contains(problem)).single;
        Flower topFlower = topFlowerVertex.wrappingFlowerList[topFlowerVertex.wrappingFlowerList.length-2];

        List<Flower> flowerPath = findEvenPath(problem.peduncle, (Flower f) => identical(f, topFlower), 1);
        Set<Dumbbell> dumbbells = new Set();

        createDumbbells() {
          Set<Flower> dumbbellFlowers = problem.flowers.difference(new Set.from(flowerPath));

          dumbbellFlowers.forEach((Flower f) => f.wrappingFlower = null);
          Set<Edge> toRemoveFromL = dumbbellFlowers.isNotEmpty ?
              dumbbellFlowers.map((Flower f) => f.edgeBorder.intersection(edgesL)).reduce((x, y) => x.union(y))
            : // The path used all of the flowers, remove just the blocking edge between the first and the last
              flowerPath.first.edgeBorder.intersection(flowerPath.last.edgeBorder).intersection(edgesL);

          while (dumbbellFlowers.isNotEmpty) {
            Edge pairingEdge = dumbbellFlowers.first.edgeBorder.intersection(edgesM).single;
            Flower f1 = pairingEdge.v1.wrappingFlowerList.last;
            Flower f2 = pairingEdge.v2.wrappingFlowerList.last;
            Dumbbell newDumbbell = new Dumbbell(f1, f2, pairingEdge);
            f1.parentStructure = newDumbbell;
            f2.parentStructure = newDumbbell;
            dumbbells.add(newDumbbell);
            assert(dumbbellFlowers.remove(f1));
            assert(dumbbellFlowers.remove(f2));
          }

          edgesL.removeAll(toRemoveFromL);

        }

        // Create HTN-s
        for (int i = 0; i < flowerPath.length; i++) {
          flowerPath[i].parentStructure = new HungarianTreeNode(flowerPath[i], null, (i + 1) % 2);
        }

        // Wire them with each other
        for (int i = flowerPath.length - 1; i > 0; i--) {
          (flowerPath[i].parentStructure as HungarianTreeNode).children.add(flowerPath[i - 1].parentStructure);
          flowerPath[i].parentStructure.children.forEach((HungarianTreeNode htn) => htn.parent = flowerPath[i].parentStructure);
        }

        printF("Flower path: $flowerPath");
        printF("Parent structure: ${flowerPath.last.parentStructure}");

        // Wire them to the tree

        // Wire the first from top - as it's on an odd level, there's always a parent
        problem.parentStructure.parent.children.remove(problem.parentStructure);
        problem.parentStructure.parent.children.add(flowerPath.last.parentStructure);
        flowerPath.last.parentStructure.parent = problem.parentStructure.parent;

        // Wire the last
        flowerPath[0].parentStructure.children.addAll(problem.parentStructure.children);
        flowerPath[0].parentStructure.children.forEach((HungarianTreeNode htn) => htn.parent = flowerPath[0].parentStructure);

        createDumbbells();

        stop(0);
        // End of P1
      } else {
        assert(problem is Edge);
        if (problem.checkP2) {
          start(1);
          // P2 happened
          printF("P2");

          HungarianTreeNode node = problem.v1.parentStructure is HungarianTreeNode ? problem.v1.parentStructure : problem.v2.parentStructure;
          Dumbbell dumbbell = problem.v1.parentStructure is Dumbbell ? problem.v1.parentStructure : problem.v2.parentStructure;

          // Part of dumbbell connected by that edge
          Flower connected = dumbbell.flowers.where((Flower f) => f.edgeBorder.contains(problem)).single;
          Flower other = dumbbell.flowers.where((Flower f) => !identical(f, connected)).single;

          // Node on odd level
          HungarianTreeNode oddNode = new HungarianTreeNode(connected, node, 1);
          connected.parentStructure = oddNode;
          // Node on even level, node from the other flower from dumbbell
          HungarianTreeNode pairNode = new HungarianTreeNode(other, oddNode, 0);
          other.parentStructure = pairNode;

          edgesL.add(problem);
          node.children.add(oddNode);
          oddNode.children.add(pairNode);
          dumbbells.remove(dumbbell);
          stop(1);
          // End of P2
        } else if (problem.checkP3) {
          start(2);
          // P3 happened
          printF("P3");

          HungarianTreeNode n1 = problem.v1.parentStructure;
          HungarianTreeNode n2 = problem.v2.parentStructure;

          List<HungarianTreeNode> pathToRoot(HungarianTreeNode htn) {
            HungarianTreeNode current = htn;
            List path = [current];
            while (!current.isRoot) {
              current = current.parent;
              path.add(current);
            }
            return path;
          }

          List path1 = new List.from(pathToRoot(n1).reversed);
          List path2 = new List.from(pathToRoot(n2).reversed);

          int i = 0;
          while ((i < path1.length) && (i < path2.length) && identical(path1[i], path2[i])) i++;

          HungarianTreeNode nearestCommonAncestor = path1[i - 1];

          Set<HungarianTreeNode> inNewFlower = new Set.from(path1.sublist(i - 1))..addAll(path2.sublist(i - 1));

          // All nodes which are not in this flower, but children to nodes in that flower
          Set<HungarianTreeNode> childrenToWire = (inNewFlower.map((HungarianTreeNode htn) => htn.children).reduce((Set x, Set y) => x.union(y)) as Set).difference(inNewFlower);

          ComplexFlower newFlower = new ComplexFlower(inNewFlower.map((HungarianTreeNode htn) => htn.node).toSet(), nearestCommonAncestor.node);

          newFlower.flowers.forEach((Flower f) {
            f.parentStructure = null;
            f.wrappingFlower = newFlower;
          });

          HungarianTreeNode newNode = new HungarianTreeNode(newFlower, nearestCommonAncestor.parent, 0);
          if (newNode.isRoot) {
            newNode.tree = nearestCommonAncestor.tree;
            newNode.tree.root = newNode;
          } else {
            nearestCommonAncestor.parent.children.remove(nearestCommonAncestor);
            nearestCommonAncestor.parent.children.add(newNode);
          }
          newFlower.parentStructure = newNode;
          newNode.children = childrenToWire;
          newNode.children.forEach((HungarianTreeNode htn) => htn.parent = newNode);

          // Debug print
          newFlower.flowers.forEach((Flower f) => printF("Parent structure of $f is ${f.parentStructure}"));

          // Add the edge to L
          edgesL.add(problem);

          stop(2);
          // End of P3
        } else if (problem.checkP4) {

          start(3);
          printF("P4");

          // P4 happened
          HungarianTreeNode n1 = problem.v1.parentStructure;
          HungarianTreeNode n2 = problem.v2.parentStructure;

          SingleVertex inN1 = problem.vertices.where((SingleVertex v) => identical(v.parentStructure, n1)).single;
          SingleVertex inN2 = problem.vertices.where((SingleVertex v) => identical(v.parentStructure, n2)).single;

          createPathNodes(HungarianTreeNode node) {
            HungarianTreeNode current = node;
            List<Flower> pathNodes = [current.node];
            while (!current.isRoot) {
              current = current.parent;
              pathNodes.add(current.node);
            }
            return pathNodes;
          }

          List<Flower> pathNodes1 = createPathNodes(n1);
          List<Flower> pathNodes2 = createPathNodes(n2);

          getPeduncleVertex(Flower f) {
            var current = f;
            while (current is ComplexFlower) current = current.peduncle;
            return current;
          }

          SingleVertex end1 = getPeduncleVertex(pathNodes1.last);
          SingleVertex end2 = getPeduncleVertex(pathNodes2.last);

          getFlowersFromEdge(Edge e) => [e.v1.wrappingFlowerList.last, e.v2.wrappingFlowerList.last];
          onThisLevel(Flower f) => f.wrappingFlowerList.last;

          Set<Edge> edges1 = makeAlterPathFromFlowers(pathNodes1, inN1, end1, getFlowersFromEdge, onThisLevel, 1).toSet();
          Set<Edge> edges2 = makeAlterPathFromFlowers(pathNodes2, inN2, end2, getFlowersFromEdge, onThisLevel, 1).toSet();

          Set<Edge> allEdges = edges1.union(edges2)..add(problem);

          // Add problem edge to L temporarily
          edgesL.add(problem);

          Set<Edge> pairingEdges = allEdges.intersection(edgesM);
          Set<Edge> blockingEdges = allEdges.intersection(edgesL);

          // Swap pairing and blocking edges
          edgesM.removeAll(pairingEdges);
          edgesL.removeAll(blockingEdges);
          edgesM.addAll(blockingEdges);
          edgesL.addAll(pairingEdges);

          // Find new peduncles in traversed flowers
          pathNodes1.forEach((Flower f) => f is ComplexFlower ? f.fixPeduncle() : null);
          pathNodes2.forEach((Flower f) => f is ComplexFlower ? f.fixPeduncle() : null);

          HungarianTree ht1 = n1.tree;
          HungarianTree ht2 = n2.tree;

          Set<Edge> allTreeEdges = ht1.collectEdges().union(ht2.collectEdges())..add(problem);
          Set<Edge> dumbbellEdges = allTreeEdges.intersection(edgesM);

          // Remove all blocking edges from the tree
          edgesL.removeAll(allTreeEdges.intersection(edgesL));

          // Create dumbbells
          dumbbellEdges.forEach((Edge e) {
            Dumbbell dumbbell = new Dumbbell(e.v1.wrappingFlowerList.last, e.v2.wrappingFlowerList.last, e);
            dumbbell.flowers.forEach((Flower f) => f.parentStructure = dumbbell);
            dumbbells.add(dumbbell);
          });

          // Trees have disappeared
          trees.remove(ht1);
          trees.remove(ht2);

          stop(3);
          // End of P4
        } else {
          // Not a problem
          printF("$problem is Not a problem anymore");

        }
      }
    }


  }

}

abstract class Flower {

  num charge = 0;

  Flower _wrappingFlower = null;

  /// Nearest enclosing flower
  Flower get wrappingFlower => _wrappingFlower;
  set wrappingFlower(Flower f) {
    _wrappingFlower = f;
    updateWrappingFlowers();
  }

  /// Set only for the outermost flower
  dynamic _parentStructure = null;

  /// Can be either an instance of HungarianTreeNode or Dumbbell
  dynamic get parentStructure => wrappingFlower == null ? _parentStructure : wrappingFlowerList.last.parentStructure;
  set parentStructure(dynamic structure) {
    wrappingFlower == null ? null : wrappingFlower = null;
    _parentStructure = structure;
  }

  Set<Edge> get edgeBorder;

  Set<Edge> get innerEdgeBorder => wrappingFlower == null ? edgeBorder : edgeBorder.difference(wrappingFlower.edgeBorder);

  bool get checkZeroCharge => true;

  _wrappingFlowerSetList(wF) {
    Flower current = this;
    while (current != null) {
      wF.add(current);
      current = current.wrappingFlower;
    }
    return wF;
  }

  _updateThisWrappingFlowers() {
    if (wrappingFlower != null) {
      _cachedWrappingFlowerList = new List.from(wrappingFlower.wrappingFlowerList)..insert(0, this);
      _cachedWrappingFlowerSet = new Set.from(wrappingFlower.wrappingFlowers)..add(this);
    } else {
      _cachedWrappingFlowerList = [this];
      _cachedWrappingFlowerSet = new Set()..add(this);
    }
  }

  updateWrappingFlowers();

  List<Flower> _cachedWrappingFlowerList = null;
  Set<Flower> _cachedWrappingFlowerSet = null;
  /// Returns a Set of flowers wrapping this flower from the innermost to the outermost (this included)
  Set<Flower> get wrappingFlowers => _cachedWrappingFlowerSet;
  /// Returns a list of flowers wrapping this flower from the innermost to the outermost (this included)
  List<Flower> get wrappingFlowerList => _cachedWrappingFlowerList;

  getMaxChargeToAdd() => 
      edgeBorder.length == 0 ? // this may happen only if the graph does not have 1-factor, and all vertices are contained in 1 flower 
          infiniteCharge()
        :
          min(edgeBorder.map((Edge e) => e.maxChargeToAdd).toList());

  getMaxChargeToSubtract();

  addCharge(num charge) => this.charge += charge;

  List<Edge> findAlterPath(SingleVertex start, SingleVertex end);

  toString();

}

class ComplexFlower extends Flower {

  Flower peduncle;

  Set<Flower> flowers;

  Set<Edge> _edgeBorder = null;

  ComplexFlower(this.flowers, this.peduncle) {
    _cachedWrappingFlowerSet = new Set()..add(this);
    _cachedWrappingFlowerList = [this];
  }

  bool get _upToDate => _edgeBorder != null;

  Set<Edge> get edgeBorder {
    if (_upToDate) {
      return _edgeBorder;
    }
    Set<Edge> border = new Set();
    // From all edges in edgeBorder-s from all flowers select only those, which are present only once
    // Those edges do not have both sides in this flower => they are in edgeBorder of this flower
    Map<Edge, num> edgeCounter = {};
    flowers.forEach((Flower f) => f.edgeBorder.forEach((Edge e) => edgeCounter[e] = (edgeCounter.containsKey(e) ? edgeCounter[e] + 1 : 1)));
    border = edgeCounter.keys.where((Edge e) => edgeCounter[e] == 1).toSet();
    _edgeBorder = border;
    return border;
  }

  getMaxChargeToSubtract() => addPotentialProblem(this, charge);

  /// After P4, peduncles of traversed flowers need to be updated
  fixPeduncle() {
    printF("Attempting to fix peduncle for ${this}, peduncle: $peduncle");
    printF("Inner blocking edges: ${getInnerBlockingEdges(peduncle)}");
    List peduncles = [];
    flowers.forEach((Flower f) {
      if (getInnerBlockingEdges(f).length == 2) {
        peduncle = f;
        peduncles.add(f);
      }
      if (f is ComplexFlower) f.fixPeduncle();
    });
    if (peduncles.length != 1) printF("WARNING! Wrong number of peduncles for ${this} ! $peduncles");
  }

  updateWrappingFlowers() {
    _updateThisWrappingFlowers();
    flowers.forEach((Flower f) => f.updateWrappingFlowers());
  }

  /// Finds an alternating path in this flower. One end has to be in a peduncle
  List<Edge> findAlterPath(SingleVertex start, SingleVertex end) {

    // Level flowers in this flower from the outermost flower
    // 0 = the flower is outer flower
    num flowerLevelFromTop = wrappingFlowerList.length - 1;
    num edgeLevelFromTop = flowerLevelFromTop + 1;

    onThisLevel(Flower lower) => lower.wrappingFlowerList[lower.wrappingFlowerList.length - 2 - flowerLevelFromTop];

    Flower endOnThisLevel = onThisLevel(end);
    Flower startOnThisLevel = onThisLevel(start);

    Flower other = null;
    SingleVertex lowPeduncle = null;
    SingleVertex lowEnd = null;
    if (identical(startOnThisLevel, peduncle)) {
      other = endOnThisLevel;
      lowPeduncle = start;
      lowEnd = end;
    }

    if (identical(endOnThisLevel, peduncle)) {
      other = startOnThisLevel;
      lowPeduncle = end;
      lowEnd = start;
    }
    // One of them has to be a peduncle
    assert(other != null);

    List<Flower> pathNodes = findEvenPath(peduncle, (Flower f) => identical(f, other), edgeLevelFromTop);

    Set<Flower> pathNodesSet = pathNodes.toSet();

    // Start / end in a flower
    SingleVertex startingVertex = lowPeduncle;
    SingleVertex endingVertex = null;

    // will first contain lists, containing either paths in flowers or the edges between those flowers
    List path = makeAlterPathFromFlowers(pathNodes, lowPeduncle, lowEnd, (Edge e) => getEdgeFlowers(e, edgeLevelFromTop), onThisLevel, 0);

    // Check if it's reversed
    if (identical(endOnThisLevel, peduncle)) path = new List.from(path.reversed);

    return path;
  }

  toString() {
    String ret = "{";
    flowers.forEach((Flower f) => ret += " ${f.toString()}${identical(peduncle, f) ? "P" : ""}");
    ret += " }";
    return ret;
  }
}

/// Creates an alternating path from given flowers [pathNodes], starting in [start] and ending in [end]
/// [getFlowersFromEdge] should specify, how to extract considered connected flowers from given edge
/// [onThisLevel] should return a flower on a considered level, containing the given vertex
/// [initialEdgeType] should be either 0 (if first edge is from L) or 1 (if first edge is from M)
List<Edge> makeAlterPathFromFlowers(List<Flower> pathNodes, SingleVertex start, SingleVertex end, List<Flower> getFlowersFromEdge(Edge e), Flower onThisLevel(SingleVertex v), num initialEdgeType) {
  SingleVertex startingVertex = start;
  SingleVertex endingVertex = null;
  List path = [];
  for (int i = 0; i < pathNodes.length - 1; i++) {

    // Edge in alternating path going out from pathNodes[i] to pathNodes[i+1]
    Edge usedEdge = pathNodes[i].innerEdgeBorder.intersection((i + initialEdgeType) % 2 == 0 ? edgesL : edgesM)
        .where((Edge e) => getFlowersFromEdge(e).contains(pathNodes[i + 1])).single;

    // SingleVertex in edge in pathNodes[i] is the ending vertex for path in flower pathNodes[i]
    endingVertex = usedEdge.vertices.where((SingleVertex v) => identical(onThisLevel(v), pathNodes[i])).single;
    path.add(pathNodes[i].findAlterPath(startingVertex, endingVertex));

    // SingleVertex on the other side of the edge is the starting vertex for path in pathNodes[i+1]
    startingVertex = usedEdge.vertices.where((SingleVertex v) => !identical(v, endingVertex)).single;
    path.add([usedEdge]);
  }
  printF("Looking for path from ${startingVertex.vertex} to ${end.vertex}");
  path.add(pathNodes.last.findAlterPath(startingVertex, end));
  return concat(path);
}

class SingleVertex extends Flower {

  Graph graph;
  num vertex;

  SingleVertex(num this.vertex) {
    updateWrappingFlowers();
  }

  bool get checkZeroCharge => false;
  Set<Edge> get edgeBorder => Graph.getNeighbors(vertex);

  getMaxChargeToSubtract() => infiniteCharge();

  updateWrappingFlowers() => _updateThisWrappingFlowers();

  List<Edge> findAlterPath(SingleVertex s1, SingleVertex s2) {
    assert(identical(s1, s2) && identical(s1, this));
    return [];
  }

  toString() => "$vertex";

}

class Dumbbell {

  Flower _f1;
  Flower _f2;

  Edge fullEdge;

  Set<Flower> get flowers => new Set.from([_f1, _f2]);

  Dumbbell(Flower this._f1, Flower this._f2, Edge this.fullEdge);

  toString() => "${_f1}<=>${_f2}";

}

class HungarianTree {

  HungarianTreeNode root;

  get maxChargeToAdd => _getMaxChargeToAdd(root);

  _getMaxChargeToAdd(HungarianTreeNode node) => min(node.children.map(_getMaxChargeToAdd).toList()..add(node.maxCharge));

  addCharge(num charge) {
    _addCharge(HungarianTreeNode node, num charge) {
      node.addCharge(charge);
      node.children.forEach((HungarianTreeNode htn) => _addCharge(htn, charge));
    }
    _addCharge(root, charge);
  }

  HungarianTree(this.root);

  Set<Edge> collectEdges() => root.collectEdges();

  toString() {
    printChildren(HungarianTreeNode htn, String offset) {
      String ret = "";
      ret += offset + htn.toString();
      htn.children.forEach((HungarianTreeNode ht) => ret += "\n" + printChildren(ht, offset + " "));
      return ret;
    }

    return "\n${printChildren(root, "")}";

  }

}

class HungarianTreeNode {

  Flower node;
  Set<HungarianTreeNode> children = new Set();

  HungarianTreeNode parent;

  /// 0 or 1
  num depthParity;
  HungarianTree _tree = null;

  get isRoot => parent == null;

  get root => isRoot ? this : parent.root;

  HungarianTree get tree => isRoot ? _tree : root.tree;
  set tree(HungarianTree ht) => isRoot ? _tree = ht : throw new Exception("Cannot change tree property in non-root node");

  get maxCharge => depthParity == 0 ? maxChargeToAdd : maxChargeToSubtract;

  get maxChargeToAdd => node.getMaxChargeToAdd();
  get maxChargeToSubtract => node.getMaxChargeToSubtract();

  addCharge(num charge) => node.addCharge(depthParity == 0 ? charge : -charge);

  Set<Edge> collectEdges() => node.edgeBorder.where((Edge e) => children.contains(e.v1.parentStructure) || children.contains(e.v2.parentStructure)).toSet().union(children.isNotEmpty ? children.map((HungarianTreeNode htn) => htn.collectEdges()).reduce((x, y) => x.union(y)) : new Set());

  HungarianTreeNode(this.node, this.parent, this.depthParity);

  toString() => "<$node>${depthParity == 1 ? "-" : "+"}";

}

File output = new File("out6.txt");
var debug = false;
printF(String s) {
  if (debug) {
    var f = output.openSync(mode: FileMode.APPEND);
    f.writeStringSync(s+"\n");
    f.closeSync();
  }
}

void main() {

  Stopwatch stopWatch = new Stopwatch();

  output.existsSync() ? output.deleteSync() : null;

  stopWatch.start();
  File file = new File("1000/1000.in");

  List<List> input = (file.readAsLinesSync()..removeWhere((String s) => s.isEmpty)).map((String s) => s.split(" ").map(int.parse).toList()).toList();

  num numVertices = input.removeAt(0)[0];

  Graph graph = new Graph(numVertices, input);
  num iteration = 0;
  try {
    while (!graph.isPerfectPairing && !graph.no1Factor) {
      iteration++;
      print("\nStarting iteration $iteration");
      print("Number of trees: ${graph.trees.length}, number of Dumbbells ${graph.dumbbells.length}");
      graph.nextIteration();
    }
    print("Success !");
    print("Price: ${edgesM.map((Edge e) => e.weight).reduce((x,y) => x+y)}");
    print("Edges: $edgesM");

  } catch (e, s) {
    printF("Error: $e, $s");
    printF("Trees: \n ${graph.trees}");
    print("Dumbbells: \n${graph.dumbbells}");
    print("Unsuccessful: \nError: $e,$s");
    print("Trees: \n${graph.trees}");
    print("Dumbbells: \n${graph.dumbbells}");
  }
  if (graph.no1Factor) print("Graph has no perfect matching");
  printF("Price: ${edgesM.map((Edge e) => e.weight).reduce((x,y) => x+y)}");
  printF("Edges: $edgesM");
  Set usedVertices = new Set();
  edgesM.forEach((Edge e) => usedVertices.addAll(e.vertices.map((SingleVertex sv) => sv.vertex)));
  if (usedVertices.length != numVertices) print("Warning: M is NOT 1 factor!");
  print("Took: ${stopWatch.elapsed}");
  print("Solutions: \n${solutionTime.map((List l) => "time: ${l[0].elapsed}, cases: ${l[1]}").reduce((x,y) => x+"\n"+y)}");
  print("Finding eps: \n${epsTime.elapsed}, cases: ${iteration}");
}
