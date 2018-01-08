package tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class Tree {

	public TreeNode root;

	private double matrix[][] = {{1.0, 0.6269006474780325, 0.8834172136770733, 0.8813063680908603, 0.8967492986790737, 0.658935434073258, 0.9550270635356867, 0.9434830515856518, 0.5642905318337071, 0.5406753639834739}, 
	{0.6269006474780325, 1.0, 0.23316997613327772, 0.7881421309498965, 0.7031296792190772, 0.24664079168440023, 0.6993937729892425, 0.978016231286417, 0.05121109036015914, 0.580515508564963},
	 {0.8834172136770733, 0.23316997613327772, 1.0, 0.38175876860601365, 0.27962882356088314, 0.4941502606295203, 0.6064421452384877, 0.694962304500226, 0.5326339016114388, 0.7928080492396583},
	 {0.8813063680908603, 0.7881421309498965, 0.38175876860601365, 1.0, 0.07192258039731547, 0.4072684451893148, 0.9552947086849037, 0.621584607555553, 0.8981593068912009, 0.33766738348691316},
	 {0.8967492986790737, 0.7031296792190772, 0.27962882356088314, 0.07192258039731547, 1.0, 0.5112541479713096, 0.9490791324473662, 0.07636916921507142, 0.9696957161299907, 0.17856974996652064},
	 {0.658935434073258, 0.24664079168440023, 0.4941502606295203, 0.4072684451893148, 0.5112541479713096, 1.0, 0.953895758657087, 0.7275465889811945, 0.6576079614232362, 0.8340566137901534}, 
	{0.9550270635356867, 0.6993937729892425, 0.6064421452384877, 0.9552947086849037, 0.9490791324473662, 0.953895758657087, 1.0, 0.7357572838276935, 0.6583479355153539, 0.9960964616492712}, 
	{0.9434830515856518, 0.978016231286417, 0.694962304500226, 0.621584607555553, 0.07636916921507142, 0.7275465889811945, 0.7357572838276935, 1.0, 0.994301338125384, 0.055722344357106435},
	 {0.5642905318337071, 0.05121109036015914, 0.5326339016114388, 0.8981593068912009, 0.9696957161299907, 0.6576079614232362, 0.6583479355153539, 0.994301338125384, 1.0, 0.47721498318270494},
	 {0.5406753639834739, 0.580515508564963, 0.7928080492396583, 0.33766738348691316, 0.17856974996652064, 0.8340566137901534, 0.9960964616492712, 0.055722344357106435, 0.47721498318270494, 1.0}};



	//private double matrix[][];

	private int size;

	/**
	 * The cluster number for each node.
	 */
	public static int[] clusterNumbers;

	public int clusterNum = 0;

	/**
	 * Is the respective node visited?
	 */
	boolean[] visited;

	public Tree(int size) {
		// this.matrix = new double[size][size];
		this.size = size;
		this.root = new TreeNode(-1);
	}

	public void buildMatrix(double[][] paraMatrix) {

		// for (int i = 0; i < size; i++) {
		// matrix[i][i] = 1;
		// for (int j = i + 1; j < size; j++) {
		// matrix[i][j] = Math.random();
		// matrix[j][i] = matrix[i][j]; // must
		// }
		//
		// }
		matrix = paraMatrix;
		//System.out.println("the ma" + Arrays.deepToString(matrix));

	}// buildMatrix

	public void buildMatrix() {

		for (int i = 0; i < size; i++) {
			matrix[i][i] = 1;
			for (int j = i + 1; j < size; j++) {
				matrix[i][j] = Math.random();
				matrix[j][i] = matrix[i][j]; // must
			}
		}
		System.out.println("the ma" + Arrays.deepToString(matrix));

	}// buildMatrix

	public void buildTree() {
		List<Integer> collections = new LinkedList<>();
		for (int i = 0; i < size; i++) {
			collections.add(i);
		}
		dfs(root,collections);
	}

	//递归划分集合建树
	private void dfs(TreeNode node,List<Integer>collections){

		if (node==null|| collections==null || collections.isEmpty())
			return;
		Node minNode=findMin(collections);

		TreeNode childNode1 = null, childNode2 = null;
		if (minNode.getX() != null) {
			childNode1 = new TreeNode(minNode.getX());
		}
		if (minNode.getY() != null) {
			childNode2 = new TreeNode(minNode.getY());
		}
		if (childNode1 == null && childNode2 == null)
			return;
		node.setChild(childNode1);
		node.setChild(childNode2);

		List<Integer> listA = new LinkedList<>();
		List<Integer> listB = new LinkedList<>();

		if (!collections.isEmpty()){
			// 划分A,B集合
			for (int i=0;i<collections.size();i++){
				if (matrix[collections.get(i)][minNode.getX()]>matrix[collections.get(i)][minNode.getY()]){
					listA.add(collections.get(i));
				}else {
					listB.add(collections.get(i));
				}
			}//of for
			dfs(childNode1,listA);
			dfs(childNode2,listB);
		}//of if
	}//of dfs

	// 层次遍历
	public void printTree() {
		if (root == null) {
			return;
		}
		ArrayList<Integer> result = new ArrayList<Integer>();
		LinkedList<TreeNode> queue = new LinkedList<TreeNode>();

		queue.offer(root);
		// result.add(root.getValue());
		while (!queue.isEmpty()) {
			int size = queue.size();
			for (int i = 0; i < size; i++) {
				TreeNode head = queue.poll();

				if (head.getLeftChild() != null) {
					queue.offer(head.getLeftChild());
					result.add(head.getLeftChild().getValue());
				}
				if (head.getRightChild() != null) {
					queue.offer(head.getRightChild());
					result.add(head.getRightChild().getValue());
				}
			}

		} // of while
		System.out.println("result" + result);
		// return result;
	} //

	public void clusteringTree(TreeNode node, List<Integer> list, int clusteringNum,double th) {
		// 最开始进来的是虚拟节点的左孩子

		//System.out.println("当前节点index" + node.getValue() + ",当前簇数" + clusteringNum);
		clusterNumbers[node.getValue()] = clusteringNum;
		list.add(node.getValue()); // 先将根节点存入list
		//System.out.println("the list" + list);
		// 如果左子树不为空继续往左找，在递归调用方法的时候一直会将子树的根存入list，这就做到了先遍历根节点,
		if ((node.getLeftChild() != null)) {
			if (matrix[node.getValue()][node.getLeftChild().getValue()] <= th) {
				clusterNum++;
				clusteringTree(node.getLeftChild(), list, clusterNum,th);
			} else {
				clusteringTree(node.getLeftChild(), list, clusteringNum,th);
			}
		}
		// 无论走到哪一层，只要当前节点左子树为空，那么就可以在右子树上遍历，保证了根左右的遍历顺序
		if (node.getRightChild() != null) {
			if (matrix[node.getValue()][node.getRightChild().getValue()] <= th) {
				clusterNum++;
				clusteringTree(node.getRightChild(), list, clusterNum,th);
			} else {
				clusteringTree(node.getRightChild(), list, clusteringNum,th);
			}

		}

	}// clusteringTree

	// 查找相似度最小的一组对象
	public Node findMin(List<Integer> list) {
		int x1 = -1, x2 = -1;
		double min = 2;

		for (int i = 0; i < list.size(); i++) {
			for (int j = i + 1; j < list.size(); j++) {
				if (matrix[list.get(i)][list.get(j)] < min) {
					min = matrix[list.get(i)][list.get(j)];
					x1 = list.get(i);
					x2 = list.get(j);
				} // of if
			} // of for j
		} // of for i

		Node node = new Node();
		if (x1 == -1) {

			node.setX(list.get(0));
			list.remove(0);
			return node;
		}

		list.remove((Object) x1);
		list.remove((Object) x2);

		node.setX(x1);
		node.setY(x2);
		return node;

	}

	public static void main(String[] args) {

		double[][] matrix = new double[10][10];

		Tree similarity = new Tree(10);
		//similarity.buildMatrix();
		similarity.buildTree();
		// similarity.printTree();
		ArrayList<Integer> list = new ArrayList<Integer>();
		clusterNumbers = new int[similarity.matrix.length];

		for (int i = 0; i < clusterNumbers.length; i++) {
			clusterNumbers[i] = -1;
		}
		int tempClusterNum = 0;

		// 遍历左树
		similarity.clusteringTree(similarity.root.getLeftChild(), list, tempClusterNum,0.6);
		System.out.println("clusterNum" + tempClusterNum);
		// 遍历右树
		// 找到当前最大簇号
		int max = Integer.MIN_VALUE;
		for (int i = 0; i < clusterNumbers.length; i++) {

			if (clusterNumbers[i] > max)
				max = clusterNumbers[i];
		} // of for i
		

		similarity.clusteringTree(similarity.root.getRightChild(), list, ++similarity.clusterNum,0.6);
		System.out.println("the clustering" + Arrays.toString(clusterNumbers));
	}

}
