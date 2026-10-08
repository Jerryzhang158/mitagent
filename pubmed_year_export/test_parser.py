import gzip
import json
import tempfile
import unittest
from pathlib import Path

from pubmed_year_export.download_pubmed_year import parse_xml_gz


SAMPLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation Status="MEDLINE">
      <PMID Version="1">12345678</PMID>
      <DateRevised><Year>2025</Year><Month>12</Month><Day>01</Day></DateRevised>
      <Article>
        <Journal>
          <ISSN IssnType="Electronic">1234-5678</ISSN>
          <JournalIssue><Volume>12</Volume><Issue>3</Issue><PubDate><Year>2025</Year></PubDate></JournalIssue>
          <Title>Example Journal</Title><ISOAbbreviation>Ex J</ISOAbbreviation>
        </Journal>
        <ArticleTitle>Effect of <i>miRNA</i> on a target.</ArticleTitle>
        <Pagination><MedlinePgn>10-20</MedlinePgn></Pagination>
        <Abstract>
          <AbstractText Label="BACKGROUND" NlmCategory="BACKGROUND">First <b>section</b>.</AbstractText>
          <AbstractText Label="RESULTS" NlmCategory="RESULTS">Second section.</AbstractText>
          <CopyrightInformation>Copyright 2025.</CopyrightInformation>
        </Abstract>
        <AuthorList><Author ValidYN="Y"><LastName>Doe</LastName><ForeName>Jane</ForeName><AffiliationInfo><Affiliation>Example University</Affiliation></AffiliationInfo></Author></AuthorList>
        <Language>eng</Language>
        <PublicationTypeList><PublicationType UI="D016428">Journal Article</PublicationType></PublicationTypeList>
      </Article>
      <OtherAbstract Type="plain-language-summary" Language="eng"><AbstractText>Plain summary.</AbstractText></OtherAbstract>
      <MeshHeadingList><MeshHeading><DescriptorName UI="D000001" MajorTopicYN="Y">Example</DescriptorName></MeshHeading></MeshHeadingList>
      <KeywordList Owner="NOTNLM"><Keyword MajorTopicYN="N">miRNA</Keyword></KeywordList>
    </MedlineCitation>
    <PubmedData><ArticleIdList><ArticleId IdType="pubmed">12345678</ArticleId><ArticleId IdType="doi">10.1/example</ArticleId><ArticleId IdType="pmc">PMC123</ArticleId></ArticleIdList></PubmedData>
  </PubmedArticle>
</PubmedArticleSet>
"""


class ParserTest(unittest.TestCase):
    def test_structured_abstract_and_mixed_content(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.xml.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(SAMPLE)
            records = parse_xml_gz(path)

        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record["pmid"], "12345678")
        self.assertEqual(record["title"], "Effect of miRNA on a target.")
        self.assertEqual(
            record["abstract"],
            "BACKGROUND: First section.\nRESULTS: Second section.",
        )
        self.assertEqual(record["abstract_sections"][0]["label"], "BACKGROUND")
        self.assertEqual(record["other_abstracts"][0]["type"], "plain-language-summary")
        self.assertEqual(record["identifiers"]["doi"], "10.1/example")
        self.assertEqual(record["mesh_headings"][0]["descriptor"], "Example")
        self.assertEqual(record["authors"][0]["affiliations"], ["Example University"])


if __name__ == "__main__":
    unittest.main()
