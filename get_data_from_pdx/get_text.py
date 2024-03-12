import pdfplumber

def extract_text_from_pdf(pdf_path):
    text = ''
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            text += page.extract_text() + '\n'
    return text

if __name__ == "__main__":
    pdf_path = 'si.pdf'
    extracted_text = extract_text_from_pdf(pdf_path)
    print(extracted_text)

