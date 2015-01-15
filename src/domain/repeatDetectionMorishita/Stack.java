package domain.repeatDetectionMorishita;

import java.util.ArrayList;

public class Stack {
	private ArrayList<ThinnerStackRecord> a;

    public Stack()
    {
       a = new ArrayList<ThinnerStackRecord>(250);       //initial capacity of 10
    }

    public boolean isEmpty()
    {
         return a.isEmpty();
    }

    public ThinnerStackRecord pop()          //pop integer element
    {
        ThinnerStackRecord last;
        last = a.remove((a.size()- 1));
        return(last);      //underflow will crash
    }

    public void push(ThinnerStackRecord x)      //push integer element
    {
        a.add(x);
    }

  
}
