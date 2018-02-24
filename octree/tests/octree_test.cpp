//
//  main.cpp
//  OctreeUnitTests
//
//  Created by Chris Satterwhite on 4/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

//#include "OctreeTests.h"
//#include "NodeTests.h"

#include <array>

#include "node.h"
#include "member.h"

#include "gtest/gtest.h"
namespace
{
  class OctreeNodeTest : public ::testing::Test
  {
    public:
      using point_type=Eigen::Vector3d;

    protected:
      OctreeNodeTest()
        : testNode(nullptr, point_type(1, 1, 1), point_type(1, 1, 1), 0, 10)
      {
      }

      node<point_type> testNode;
  };

  // Test  a simple construction
  TEST_F(OctreeNodeTest, Construction)
  {

    // Test origin is set correctly
    bool flag = true;
    for (int i=0; i<3; i++)
    {
      // Test origin is set correctly
      EXPECT_EQ(testNode.getOrigin()[i], 1);
      // Test halfDimension is set correctly
      EXPECT_EQ(testNode.getHalfDimension()[i], 1);
    }

    // Test children are NULL
    flag = true;
    for (int i=0; i<8; i++)
    {
      EXPECT_EQ(testNode.getChild(i), nullptr);
    }

    // Test level set
    EXPECT_EQ(testNode.getLevel(), 1);
  }

  // Test the leaf node feature
  TEST_F(OctreeNodeTest, LeafNode)
  {
    EXPECT_TRUE(testNode.isLeafNode());
  }

  // Test adding a member
  TEST_F(OctreeNodeTest, AddMember)
  {
    int maxMembers = 10;
    testNode.setMaxMembers(maxMembers);

    int nX = 3;
    int nY = 3;
    int nZ = 3;
    int counter = 0;
    for (int i=0; i<nX; i++)
    {
      for (int j=0; j<nY; j++)
      {
        for (int k=0; k<nZ; k++)
        {
          point_type obj; obj << i, j, k;
          member<point_type> member(&obj,obj);
          testNode.addMember(member);
          counter++;
          if (counter <= maxMembers)
          {
              EXPECT_TRUE(testNode.isLeafNode());
          }
          else
          {
              // Once maxMembers is exceeded, all members should be pushed to
              // children and therefore node will no longer be leaf and should
              // contain no members.
              EXPECT_FALSE(testNode.isLeafNode());
          }
          // Number of members in node should be the same as the number added
          // due to recursive getMembers method.
          EXPECT_EQ(testNode.getMembers().size(), counter);
        }
      }
    }
  }

  // Test finding node containing specified member
  TEST_F(OctreeNodeTest, NodeContainingMember)
  {
    using member_type=member<point_type>;

    int maxMembers = 10;
    testNode.setMaxMembers(maxMembers);

    int nX = 3;
    int nY = 3;
    int nZ = 3;
    int counter = 0;
    for (int i=0; i<nX; i++)
    {
      for (int j=0; j<nY; j++)
      {
        for (int k=0; k<nZ; k++)
        {
          point_type obj; obj << i, j, k;
          member<point_type> member(&obj,obj);
          testNode.addMember(member);
          counter++;
        }
      }
    }

    const point_type origin = testNode.getOrigin();
    point_type obj;
    obj[0] = origin[0]-1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]-1;
    member_type member0(&obj,obj);
    testNode.addMember(member0);
    EXPECT_EQ(testNode.getChildContainingMember(member0), 0);

    obj[0] = origin[0]-1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]+1;
    member_type member1(&obj,obj);
    testNode.addMember(member1);
    EXPECT_EQ(testNode.getChildContainingMember(member1), 1);

    obj[0] = origin[0]-1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]-1;
    member_type member2(&obj,obj);
    testNode.addMember(member2);
    EXPECT_EQ(testNode.getChildContainingMember(member2), 2);

    obj[0] = origin[0]-1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]+1;
    member_type member3(&obj,obj);
    testNode.addMember(member3);
    EXPECT_EQ(testNode.getChildContainingMember(member3), 3);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]-1;
    member_type member4(&obj,obj);
    testNode.addMember(member4);
    EXPECT_EQ(testNode.getChildContainingMember(member4), 4);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]+1;
    member_type member5(&obj,obj);
    testNode.addMember(member5);
    EXPECT_EQ(testNode.getChildContainingMember(member5), 5);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]-1;
    member_type member6(&obj,obj);
    testNode.addMember(member6);
    EXPECT_EQ(testNode.getChildContainingMember(member6), 6);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]+1;
    member_type member7(&obj,obj);
    testNode.addMember(member7);
    EXPECT_EQ(testNode.getChildContainingMember(member7), 7);
  }

  TEST_F(OctreeNodeTest, CreatePoint)
  {
    using member_type=member<point_type>;

    using member_type=member<point_type>;

    int maxMembers = 10;
    testNode.setMaxMembers(maxMembers);

    int nX = 3;
    int nY = 3;
    int nZ = 3;
    int counter = 0;
    for (int i=0; i<nX; i++)
    {
      for (int j=0; j<nY; j++)
      {
        for (int k=0; k<nZ; k++)
        {
          point_type obj; obj << i, j, k;
          member<point_type> member(&obj,obj);
          testNode.addMember(member);
          counter++;
        }
      }
    }

    const point_type origin = testNode.getOrigin();
    point_type obj;
    obj[0] = origin[0]-1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]-1;
    member_type member0(&obj,obj);
    testNode.addMember(member0);

    obj[0] = origin[0]-1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]+1;
    member_type member1(&obj,obj);
    testNode.addMember(member1);

    obj[0] = origin[0]-1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]-1;
    member_type member2(&obj,obj);
    testNode.addMember(member2);

    obj[0] = origin[0]-1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]+1;
    member_type member3(&obj,obj);
    testNode.addMember(member3);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]-1;
    member_type member4(&obj,obj);
    testNode.addMember(member4);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]-1;
    obj[2] = origin[2]+1;
    member_type member5(&obj,obj);
    testNode.addMember(member5);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]-1;
    member_type member6(&obj,obj);
    testNode.addMember(member6);

    obj[0] = origin[0]+1;
    obj[1] = origin[1]+1;
    obj[2] = origin[2]+1;
    member_type member7(&obj,obj);
    testNode.addMember(member7);

    testNode.createParent(0);
    EXPECT_FALSE(testNode.getChild(7) == nullptr);
    EXPECT_EQ(testNode.getLevel(), 2);
  }
}

#if 0

int main()
{
    Test::Suite tests;
    tests.add(std::auto_ptr<Test::Suite>(new NodeTests));

    Test::TextOutput output(Test::TextOutput::Verbose);
    return tests.run(output);
#endif // 0

