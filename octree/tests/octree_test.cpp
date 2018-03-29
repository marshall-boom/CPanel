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
#include <memory>

#include "node.h"
#include "member.h"
#include "octree.h"

#include "gtest/gtest.h"

namespace
{
  class OctreeNodeTest : public ::testing::Test
  {
    public:
      using point_type=Eigen::Vector3d;

    protected:
      OctreeNodeTest()
        : testNode(nullptr, point_type(1, 1, 1), point_type(1, 1, 1), 0, 10, 0.1)
      {
      }

    protected:
      using node_type = node<point_type>;
      using node_index_type = node_type::member_index_type;
      node<point_type> testNode;
  };

  // Test  a simple construction
  TEST_F(OctreeNodeTest, Construction)
  {

    // Test origin is set correctly
    for (int i=0; i<3; i++)
    {
      // Test origin is set correctly
      EXPECT_EQ(testNode.getOrigin()[i], 1);
      // Test halfDimension is set correctly
      EXPECT_EQ(testNode.getHalfDimension()[i], 1);
    }

    // Test children are NULL
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
    node_index_type maxMembers = 10;
    testNode.setMaxMembers(maxMembers);

    int nX = 3;
    int nY = 3;
    int nZ = 3;
    node_index_type counter = 0;
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
    node_index_type counter = 0;
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
    node_index_type counter = 0;
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

  class test_object
  {
    public:
      test_object(Eigen::Vector3d point) : center(point) {}

      Eigen::Vector3d getCenter() const {return center;}

    private:
      Eigen::Vector3d center;
  };

  class test_octree_class : public octree<test_object>
  {
    public:
      Eigen::Vector3d findRefPoint(const test_object &member) {return member.getCenter();}

      test_octree_class() : octree() {}
  };

  TEST(OctreeTest, Construction)
  {
    test_octree_class testOctree;
    EXPECT_EQ(testOctree.getMaxMembersPerNode(), 10);
  }

  TEST(OctreeTest, MaxMembers)
  {
    int max = 5;
    test_octree_class testOctree;

    testOctree.setMaxMembers(max);
    EXPECT_EQ(testOctree.getMaxMembersPerNode(), max);
  }

  TEST(OctreeTest, AddData)
  {
    test_octree_class testOctree;
    int nX = 3;
    int nY = 3;
    int nZ = 3;

    std::vector<test_object *> data;

    for (int i=0; i<nX; i++)
    {
      for (int j=0; j<nY; j++)
      {
        for (int k=0; k<nZ; k++)
        {
          Eigen::Vector3d center;
          center[0] = i;
          center[1] = j;
          center[2] = k;
          data.push_back(new test_object(center));
        }
      }
    }
    testOctree.addData(data);
    for (int i=0; i<nX*nY*nZ; ++i)
    {
      delete data[i];
      data[i] = nullptr;
    }

    EXPECT_EQ(testOctree.getMembers().size(), static_cast<test_octree_class::member_index_type>(nX*nY*nZ));

    for (int i=0; i<3; i++)
    {
      EXPECT_EQ(testOctree.getRootNode()->getOrigin()[i], 1);
    }
    EXPECT_TRUE(testOctree.getRootNode()->getParent() == NULL);

    node<test_object> * oldRoot = testOctree.getRootNode();

    test_object *newData;
    Eigen::Vector3d newPoint; newPoint << -2,-2,-2;
    newData = new test_object(newPoint);
    testOctree.addData(newData);

    EXPECT_EQ(testOctree.getMembers().size(), static_cast<test_octree_class::member_index_type>(nX*nY*nZ+1));
    EXPECT_FALSE(testOctree.getRootNode() == oldRoot);
  }

}
#if 0

void OctreeTests::test_addData()
{
    int nX = 3;
    int nY = 3;
    int nZ = 3;
    testObj* obj;
    for (int i=0; i<nX; i++)
    {
        for (int j=0; j<nY; j++)
        {
            for (int k=0; k<nZ; k++)
            {
                std::array<double,3> center;
                center[0] = i;
                center[1] = j;
                center[2] = k;
                obj = new testObj(center);
                data.push_back(obj);
            }
        }
    }
    testOctreeClass testOctree;
    testOctree.addData(data);

    TEST_ASSERT(testOctree.getMembers().size()==nX*nY*nZ);

    bool flag = true;
    for (int i=0; i<3; i++)
    {
        if (testOctree.getRootNode()->getOrigin()[i] != 1)
        {
            flag = false;
        }
    }
    TEST_ASSERT(flag)
    TEST_ASSERT(testOctree.getRootNode()->getParent() == NULL);

    node<testObj>* oldRoot = testOctree.getRootNode();

    testObj* newData;
    std::array<double,3> newPoint = {-2,-2,-2};
    newData = new testObj(newPoint);
    data.push_back(newData);
    testOctree.addData(newData);

    TEST_ASSERT(testOctree.getMembers().size()==nX*nY*nZ+1);
    TEST_ASSERT(testOctree.getRootNode() != oldRoot);
}
#endif // 0
