package algorithms.curves;

public class FailedToConvergeException extends Exception {

    protected static final long serialVersionUID = 456789123;

    public FailedToConvergeException(String msg) {
        super(msg);
    }
}
