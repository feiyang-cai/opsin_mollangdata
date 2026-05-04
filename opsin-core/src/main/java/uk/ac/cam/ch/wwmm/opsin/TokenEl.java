package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

class TokenEl extends Element {
	
	private String value;
	private Fragment frag;
	// MolLangData: Support for child TokenEl elements
	private final List<Element> children = new ArrayList<>();

	TokenEl(String name) {
		super(name);
		this.value = "";
	}
	
	TokenEl(String name, String value) {
		super(name);
		this.value = value;
	}

	@Override
	void addChild(Element child) {
		// MolLangData: Allow TokenEl children
		if (!(child instanceof TokenEl)) {
			throw new UnsupportedOperationException("TokenEl can only have TokenEl children");
		}
		child.setParent(this);
		children.add(child);
	}
	
	@Override
	Element copy() {
		TokenEl copy = new TokenEl(this.name, this.value);
		for (int i = 0, len = this.attributes.size(); i < len; i++) {
			Attribute atr = this.attributes.get(i);
			copy.addAttribute(new Attribute(atr));
		}
		// MolLangData: Copy children
		for (Element childEl : this.children) {
			Element newChild = childEl.copy();
			newChild.setParent(copy);
			copy.addChild(newChild);
		}
		return copy;
	}
	
	/**
	 * Creates a copy with no parent
	 * The provided value is used instead of the Element to be copied's value
	 * @param value
	 * @return
	 */
	TokenEl copy(String value) {
		TokenEl copy = new TokenEl(this.name, value);
		for (int i = 0, len = this.attributes.size(); i < len; i++) {
			Attribute atr = this.attributes.get(i);
			copy.addAttribute(new Attribute(atr));
		}
		// MolLangData: Copy children
		for (Element childEl : this.children) {
			Element newChild = childEl.copy();
			newChild.setParent(copy);
			copy.addChild(newChild);
		}
		return copy;
	}
	
	@Override
	Element getChild(int index) {
		return children.get(index);
	}

	@Override
	int getChildCount() {
		return children.size();
	}

	@Override
	List<Element> getChildElements() {
		return new ArrayList<>(children);
	}

	@Override
	List<Element> getChildElements(String name) {
		List<Element> elements = new ArrayList<>(1);
		for (Element element : children) {
			if (element.name.equals(name)) {
				elements.add(element);
			}
		}
		return elements;
	}

	@Override
	Element getFirstChildElement(String name) {
		for (Element child : children) {
			if (child.getName().equals(name)) {
				return child;
			}
		}
		return null;
	}
	
	@Override
	Element getLastChildElement() {
		int childCount = children.size();
		return childCount > 0 ? children.get(childCount - 1) : null;
	}
	
	@Override
	Fragment getFrag() {
		return frag;
	}
	
	String getValue() {
		return value;
	}

	@Override
	int indexOf(Element child) {
		return children.indexOf(child);
	}

	@Override
	void insertChild(Element child, int index) {
		// MolLangData: Allow TokenEl children
		if (!(child instanceof TokenEl)) {
			throw new UnsupportedOperationException("TokenEl can only have TokenEl children");
		}
		child.setParent(this);
		children.add(index, child);
	}
	
	@Override
	boolean removeChild(Element child) {
		child.setParent(null);
		return children.remove(child);
	}
	
	@Override
	Element removeChild(int index) {
		Element removed = children.remove(index);
		removed.setParent(null);
		return removed;
	}

	@Override
	void replaceChild(Element oldChild, Element newChild) {
		// MolLangData: Allow TokenEl children
		if (!(newChild instanceof TokenEl)) {
			throw new UnsupportedOperationException("TokenEl can only have TokenEl children");
		}
		int index = indexOf(oldChild);
		if (index == -1) {
			throw new RuntimeException("oldChild is not a child of this element.");
		}
		removeChild(index);
		insertChild(newChild, index);
	}
	
	@Override
	void setFrag(Fragment frag) {
		this.frag = frag;
	}

	void setValue(String text) {
		this.value = text;
	}

}
