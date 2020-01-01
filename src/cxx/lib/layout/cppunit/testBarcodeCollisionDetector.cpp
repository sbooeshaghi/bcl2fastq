/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file testBarcodeCollisionDetector.cpp
 *
 * \brief BarcodeCollisionDetector cppunit test declarations.
 *
 * \author Aaron Day
 */


#include "layout/BarcodeCollisionDetector.hh"

#include "RegistryName.hh"
#include "testBarcodeCollisionDetector.hh"


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TestBarcodeCollisionDetector, registryName("BarcodeCollisionDetector"));


void TestBarcodeCollisionDetector::setUp()
{
    componentMaxMismatches_.push_back(0);
    componentMaxMismatches_.push_back(1);
    componentMaxMismatches_.push_back(2);
}

void TestBarcodeCollisionDetector::tearDown()
{
}

void TestBarcodeCollisionDetector::testSuccess()
{
    bcl2fastq::layout::BarcodesContainer barcodeContainer1;
    std::vector<std::string> barcodeComponents1;
    barcodeComponents1.push_back("ACGTACGT");
    barcodeComponents1.push_back("TTTTCCCC");
    barcodeComponents1.push_back("CGCGCGCGCGCG");
    barcodeComponents1.push_back("AAAAAAAAAA");
    barcodeContainer1.push_back(bcl2fastq::layout::Barcode(barcodeComponents1.begin(), barcodeComponents1.end()));

    bcl2fastq::layout::BarcodesContainer barcodeContainer2;
    std::vector<std::string> barcodeComponents2;
    barcodeComponents2.push_back("ACGGACCT"); // 2 mismatches
    barcodeComponents2.push_back("ATAACCCC"); // 3 mismatches 
    barcodeComponents2.push_back("CACGAGTGCCCA"); // 5 mismatches
    barcodeComponents2.push_back("GGGGGGGGGG"); // all mismatches
    barcodeContainer2.push_back(bcl2fastq::layout::Barcode(barcodeComponents2.begin(), barcodeComponents2.end()));

    bcl2fastq::layout::BarcodeCollisionDetector collisionDetector(componentMaxMismatches_);

    collisionDetector.validateBarcode(barcodeContainer1);

    collisionDetector.validateBarcode(barcodeContainer2);
}

void TestBarcodeCollisionDetector::testMismatchLength()
{
    bcl2fastq::layout::BarcodesContainer barcodeContainer1;
    std::vector<std::string> barcodeComponents1;
    barcodeComponents1.push_back("ACGTACGT");
    barcodeComponents1.push_back("TTTTCCCCG");
    barcodeContainer1.push_back(bcl2fastq::layout::Barcode(barcodeComponents1.begin(), barcodeComponents1.end()));

    bcl2fastq::layout::BarcodesContainer barcodeContainer2;
    std::vector<std::string> barcodeComponents2;
    barcodeComponents2.push_back("CGCGTATA");
    barcodeComponents2.push_back("AACCGGTT");
    barcodeContainer2.push_back(bcl2fastq::layout::Barcode(barcodeComponents2.begin(), barcodeComponents2.end()));

    bcl2fastq::layout::BarcodeCollisionDetector collisionDetector(componentMaxMismatches_);

    collisionDetector.validateBarcode(barcodeContainer1);

    CPPUNIT_ASSERT_THROW(collisionDetector.validateBarcode(barcodeContainer2),
                         bcl2fastq::layout::BarcodeCollisionError);

    // 2 barcodes of different length in the same sample
    bcl2fastq::layout::BarcodesContainer barcodeContainer3;
    std::vector<std::string> barcodeComponents3;
    barcodeComponents3.push_back("AAAAAAAA");
    barcodeComponents3.push_back("TTTTTTTT");
    barcodeContainer3.push_back(bcl2fastq::layout::Barcode(barcodeComponents3.begin(), barcodeComponents3.end()));

    barcodeComponents3.clear();
    barcodeComponents3.push_back("AAAAAAAAAA");
    barcodeComponents3.push_back("TTTTTTTTTT");
    barcodeContainer3.push_back(bcl2fastq::layout::Barcode(barcodeComponents3.begin(), barcodeComponents3.end()));

    CPPUNIT_ASSERT_THROW(collisionDetector.validateBarcode(barcodeContainer3),
                         bcl2fastq::layout::BarcodeCollisionError);
}

void TestBarcodeCollisionDetector::testNumComponents()
{
    bcl2fastq::layout::BarcodesContainer barcodeContainer1;
    std::vector<std::string> barcodeComponents1;
    barcodeComponents1.push_back("ACGTACGT");
    barcodeComponents1.push_back("TTTTCCCC");
    barcodeComponents1.push_back("CTCTGATG");
    barcodeContainer1.push_back(bcl2fastq::layout::Barcode(barcodeComponents1.begin(), barcodeComponents1.end()));

    bcl2fastq::layout::BarcodesContainer barcodeContainer2;
    std::vector<std::string> barcodeComponents2;
    barcodeComponents2.push_back("CGCGTATA");
    barcodeComponents2.push_back("AACCGGTT");
    barcodeContainer2.push_back(bcl2fastq::layout::Barcode(barcodeComponents2.begin(), barcodeComponents2.end()));

    bcl2fastq::layout::BarcodeCollisionDetector collisionDetector(componentMaxMismatches_);

    collisionDetector.validateBarcode(barcodeContainer1);

    CPPUNIT_ASSERT_THROW(collisionDetector.validateBarcode(barcodeContainer2),
                         bcl2fastq::layout::BarcodeCollisionError);
}

void TestBarcodeCollisionDetector::testCollisionSameSample()
{
    // Add 2 close barcodes, which is not ok even though they're for the same sample
    bcl2fastq::layout::BarcodesContainer barcodeContainer3;
    std::vector<std::string> barcodeComponents3;
    barcodeComponents3.push_back("ACGTACGT"); // 0 mismatch
    barcodeComponents3.push_back("TTTTCCCC"); // 0 mismatches
    barcodeComponents3.push_back("CGCGCGCGCGCG"); // 0 mismatches
    barcodeComponents3.push_back("AAAAGGGGGG"); // 6 mismatches
    barcodeContainer3.push_back(bcl2fastq::layout::Barcode(barcodeComponents3.begin(), barcodeComponents3.end()));

    barcodeComponents3.clear();
    barcodeComponents3.push_back("ACGTACGT"); // 0 mismatches
    barcodeComponents3.push_back("TTTTCCCC"); // 0 mismatches
    barcodeComponents3.push_back("CGCGCGCGCGCG"); // 0 mismatches
    barcodeComponents3.push_back("AAAAAGGGGG"); // 5 mismatches
    barcodeContainer3.push_back(bcl2fastq::layout::Barcode(barcodeComponents3.begin(), barcodeComponents3.end()));

    bcl2fastq::layout::BarcodeCollisionDetector collisionDetector(componentMaxMismatches_);

    CPPUNIT_ASSERT_THROW(collisionDetector.validateBarcode(barcodeContainer3),
                         bcl2fastq::layout::BarcodeCollisionError);
}

void TestBarcodeCollisionDetector::testCollision()
{
    bcl2fastq::layout::BarcodesContainer barcodeContainer1;
    std::vector<std::string> barcodeComponents1;
    barcodeComponents1.push_back("ACGTACGT");
    barcodeComponents1.push_back("TTTTCCCC");
    barcodeContainer1.push_back(bcl2fastq::layout::Barcode(barcodeComponents1.begin(), barcodeComponents1.end()));

    bcl2fastq::layout::BarcodesContainer barcodeContainer2;
    std::vector<std::string> barcodeComponents2;
    barcodeComponents2.push_back("ACGTACGT"); // 2 mismatches
    barcodeComponents2.push_back("TTTTTTCC"); // 2 mismatches
    barcodeContainer2.push_back(bcl2fastq::layout::Barcode(barcodeComponents2.begin(), barcodeComponents2.end()));

    bcl2fastq::layout::BarcodeCollisionDetector collisionDetector(componentMaxMismatches_);

    collisionDetector.validateBarcode(barcodeContainer1);

    CPPUNIT_ASSERT_THROW(collisionDetector.validateBarcode(barcodeContainer2),
                         bcl2fastq::layout::BarcodeCollisionError);

    barcodeComponents1.clear();
    barcodeComponents1.push_back("GGGGGGGGGGGG");
    barcodeContainer1.push_back(bcl2fastq::layout::Barcode(barcodeComponents1.begin(), barcodeComponents1.end()));

    barcodeComponents2.clear();
    barcodeComponents2.push_back("AAAAGGGGGGGG"); // 4 mismatches
    barcodeContainer2.push_back(bcl2fastq::layout::Barcode(barcodeComponents2.begin(), barcodeComponents2.end()));
    CPPUNIT_ASSERT_THROW(collisionDetector.validateBarcode(barcodeContainer2),
                         bcl2fastq::layout::BarcodeCollisionError);
}

