package tree;

public class TreeNode {

    private  Integer value;//节点的下标
    
    private TreeNode rightChild;

    private TreeNode leftChild;


    public TreeNode(Integer value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public void setValue(int value) {
        this.value = value;
    }

    public TreeNode getRightChild() {
        return rightChild;
    }

    public void setRightChild(TreeNode rightChild) {
        this.rightChild = rightChild;
    }

    public TreeNode getLeftChild() {
        return leftChild;
    }

    public void setLeftChild(TreeNode leftChild) {
        this.leftChild = leftChild;
    }

    public  void setChild(TreeNode child){
        if (child==null)
            return;
        if (leftChild==null) {
            leftChild = child;
            return;
        }
        if (rightChild==null){
            rightChild=child;
            return;
        }

        // 抛异常提醒

    }
}
